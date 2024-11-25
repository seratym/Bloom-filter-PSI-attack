use std::collections::BTreeSet;
use std::{fmt, ops::Range, thread, time::Instant, vec};
use std::fmt::{Display, Formatter};
use bytevec::ByteEncodable;
use combinations::compute_unique;
use num_bigint::BigUint;
use rand::{seq::index::sample, thread_rng, Rng};
use structopt::StructOpt;
use xxh3::hash64_with_seed;
use rayon::prelude::*;
use cached::proc_macro::cached;
use topk::FilteredSpaceSaving;

mod combinations;

#[derive(StructOpt, Debug)]
#[structopt(
    name = "BF-MPSI attack",
    about = "An attack on Bloom filter-based (M)PSI protocols"
)]

struct Opt {
    /// log2 of the total number of elements in the universe
    #[structopt(short, long)]
    universe_size: usize,

    #[structopt(long)]
    universe_client_factor: usize,

    /// log2 of the max private set size
    #[structopt(short = "s", long)]
    max_set_size: usize,

    /// log2 of the client's input set.
    #[structopt(long)]
    client_set_size: usize,

    /// ln of the false positive probability power
    #[structopt(short = "e", long)]
    false_positives: isize,

    /// Whether the attack should be hidden for the client
    #[structopt(long)]
    covert: bool,

    /// Number of trials that experiment will be run with and without the target set
    #[structopt(short, long, default_value = "100000")]
    trials: usize,

    /// Verbosity of the program. -v for parameters, -vv for attack search, and -vvv for thread info
    #[structopt(short, parse(from_occurrences))]
    verbose: u8,

}

pub fn hash_element(element: &usize, seed: u64) -> usize {
    let element_bytes = (*element as u64).encode::<u64>().unwrap();
    hash64_with_seed(&element_bytes, seed) as usize
}

#[derive(Debug, Clone)]
struct BloomFilterParams {
    bin_count: usize,
    seeds: Vec<u64>,
}

impl BloomFilterParams {
    fn new(element_count: usize, e_pow: isize) -> Self {
        let mut rng = thread_rng();
        let hash_count = (-e_pow) as usize;

        let funny_looking_thing = 2f64.powf(e_pow as f64 / hash_count as f64);
        let bin_count = (-(hash_count as f64) * (element_count as f64 + 0.5) / (1. - funny_looking_thing).ln()).ceil() + 1.;
        // let seeds: Vec<u64> = (0..hash_count).map(|_| rng.gen()).collect();
        let seeds: Vec<u64> = (0..hash_count).map(|x| x as u64).collect();
        BloomFilterParams {
            bin_count: bin_count as usize,
            seeds,
        }
    }
}

fn hash_hits(bf_params: &BloomFilterParams, set: &Range<usize>) -> Vec<Vec<usize>> {
    let mut result: Vec<Vec<usize>> = (0..bf_params.bin_count).map(|_| vec![]).collect();

    for el in set.clone() {
        for seed in bf_params.seeds.iter() {
            result[hash_element(&el, *seed) % bf_params.bin_count].push(el);
        }
    }

    result
}

fn find_non_unique(bins: &[BTreeSet<usize>]) -> Option<usize> {
    let mut seen_elements = BTreeSet::new();

    for bin in bins {
        for element in bin {
            if seen_elements.contains(element) {
                return Some(*element);
            }
            seen_elements.insert(element);
        }
    }

    None
}

#[cached]
fn compute_exact_score(universe_size: usize, page_size: usize, bins: Vec<BTreeSet<usize>>) -> BigUint {
    if universe_size == 0 {
        return BigUint::ZERO;
    }

    // Check if any bin is empty
    for bin in &bins {
        if bin.is_empty() {
            return BigUint::ZERO;
        }
    }

    let split_element = find_non_unique(&bins);

    match split_element {
        Some(element) => {
            // There is non-unique element
            let mut bins_remove_element = vec![];
            let mut bins_pick_element = vec![];
            for bin in bins {
                // Remove
                let mut new_bin = bin.clone();
                new_bin.remove(&element);
                bins_remove_element.push(new_bin);

                // Pick
                if !bin.contains(&element) {
                    bins_pick_element.push(bin);
                }
            }

            // Compute total
            let score_remove = compute_exact_score(universe_size - 1, page_size, bins_remove_element);
            let score_pick = compute_exact_score(universe_size - 1, page_size - 1, bins_pick_element);

            score_remove + score_pick
        },
        None => {
            // All bins contain unique elements
            let bin_sizes: Vec<usize> = bins.iter().map(|bin| bin.len()).collect();
            compute_unique(&bin_sizes, page_size, universe_size)
        }
    }
}

fn compute_approx_score_log2(universe_size: usize, page_size: usize, bins: Vec<BTreeSet<usize>>) -> u64 {
    compute_exact_score(universe_size, page_size, bins).bits()
}

fn snipe(bf_params: &BloomFilterParams, u1: Range<usize>, x1_size: usize, u2: Range<usize>, x2_size: usize, thread_count: usize, covert: bool, verbose: u8) -> (Vec<usize>, usize) {
    let hh = hash_hits(bf_params, &u2);

    // Find target
    let (target, avoid_bin) = find_target(bf_params, &u2, &hh);
    if verbose >= 1 { println!("% Targetting element {:?}", target) };

    // Find bin to snipe the target with
    let bins_target: Vec<usize> = bf_params.seeds.iter()
        .map(|seed| hash_element(&target, *seed) % bf_params.bin_count)
        .collect();
    let bins_target: BTreeSet<usize> = BTreeSet::from_iter(bins_target);

    let u2_len = u2.len();

    // First filter away any element in the universe that do not snipe the target
    let all_candidates: Vec<usize> = find_all_candidates(bf_params, u1, thread_count, covert, &hh, avoid_bin, &bins_target);

    if verbose >= 1 { println!("% Found {:?} candidates", all_candidates.len()) };

    // Create attack set
    if all_candidates.len() <= x1_size {
        return (all_candidates, target)
    }

    let x1 = search_best_candidates(bf_params, all_candidates, x1_size, u2_len, x2_size, &hh, thread_count);
    (x1, target)
}

fn search_best_candidates(bf_params: &BloomFilterParams, all_candidates: Vec<usize>, x1_size: usize, u2_len: usize, x2_size: usize, hh: &Vec<Vec<usize>>, thread_count: usize) -> Vec<usize> {
    let mut topk = FilteredSpaceSaving::new(x1_size);
    let chunk_size = (all_candidates.len() + thread_count - 1)/thread_count; // ceil the chunk_size

    thread::scope(|s| {
        let mut handles = Vec::with_capacity(thread_count);

        let hh_ref = &hh;
        for candidate_chunk in all_candidates.chunks(chunk_size) {
            handles.push(s.spawn(move || {
                let mut topk: FilteredSpaceSaving<usize> = FilteredSpaceSaving::new(x1_size);
                for &candidate in candidate_chunk {
                    let candidate_indices: BTreeSet<usize> = BTreeSet::from_iter(bf_params.seeds.iter().map(|&seed| hash_element(&candidate, seed) % bf_params.bin_count));
                    let bins: Vec<BTreeSet<usize>> = candidate_indices.iter().map(|&index|
                        BTreeSet::from_iter(hh_ref[index].iter().copied())
                    ).collect();
                    let score = compute_approx_score_log2(u2_len, x2_size, bins);
                    topk.insert(candidate, score);
                }

                topk
            }));
        }

        for handle in handles {
            let sub_topk = handle.join().expect("Failed to join thread");
            topk.merge(&sub_topk).expect("Merging errored");
        }
    });
    topk.iter().map(|(candidate, _)| *candidate).collect()
}

fn find_all_candidates(bf_params: &BloomFilterParams, u1: Range<usize>, thread_count: usize, covert: bool, hh: &Vec<Vec<usize>>, avoid_bin: usize, bins_target: &BTreeSet<usize>) -> Vec<usize> {
    let mut all_candidates = vec![];
    let chunk_size = u1.len() / thread_count;

    thread::scope(|s| {
        let mut handles = Vec::with_capacity(thread_count);

        let hh_ref = &hh;
        let bins_target_ref = &bins_target;

        for i in 0..thread_count {
            handles.push(s.spawn(move || {
                let mut candidates = vec![];
                for y in (u1.start + i * chunk_size)..((i + 1) * chunk_size) {
                    let sniped_target;
                    if covert {
                        sniped_target = potential_candidate(bf_params, hh_ref, bins_target_ref, &y, avoid_bin);
                    } else {
                        sniped_target = potential_candidate(bf_params, hh_ref, bins_target_ref, &y, usize::MAX);
                    }
                    if sniped_target {
                        candidates.push(y);
                    }
                }

                candidates
            }));
        }

        for handle in handles {
            let candidates = handle.join().expect("Failed to join thread");
            all_candidates.extend(candidates);
        }
    });
    all_candidates
}

fn potential_candidate(bf_params: &BloomFilterParams, hh_ref: &Vec<Vec<usize>>, bins_target_ref: &BTreeSet<usize>, y: &usize, avoid_bin: usize) -> bool {
    let mut sniped_target = false;

    for seed in bf_params.seeds.iter() {
        let bloomf_index = hash_element(&y, *seed) % bf_params.bin_count;

        if hh_ref[bloomf_index].is_empty() || bloomf_index == avoid_bin {
            return false;
        } else if hh_ref[bloomf_index].len() == 1 && bins_target_ref.contains(&bloomf_index) {
            sniped_target = true;
        }
    }
    sniped_target
}

fn find_target(bf_params: &BloomFilterParams, u2: &Range<usize>, hash_hits: &Vec<Vec<usize>>) -> (usize, usize) {
    let mut best_score = usize::MAX;
    let mut best_x: usize = usize::MAX;

    for x in u2.clone() {
        let score: usize = bf_params.seeds.iter()
            .map(|seed| hash_hits[hash_element(&x, *seed) % bf_params.bin_count].len())
            .sum();

        if score < best_score {
            best_score = score;
            best_x = x;
        }
    }

    let avoid_bin = bf_params.seeds.iter()
        .map(|seed| {
            let index = hash_element(&best_x, *seed) % bf_params.bin_count;
            (index, hash_hits[index].len())
        })
        .min_by_key(|&(_, len)| {
            if len == 1 {
                usize::MAX
            } else {
                len
            }
        })
        .map(|(index, _)| index);

    (best_x, avoid_bin.unwrap())
}

struct BloomFilter {
    bins: Vec<bool>,
    seeds: Vec<u64>,
}

impl BloomFilter {
    fn new(params: &BloomFilterParams) -> Self {
        BloomFilter {
            bins: vec![false; params.bin_count],
            seeds: params.seeds.clone(),
        }
    }

    fn insert(&mut self, element: usize) {
        let bin_count = self.bins.len();
        for seed in self.seeds.iter() {
            self.bins[hash_element(&element, *seed) % bin_count] = true;
        }
    }

    fn clear(&mut self) {
        self.bins = vec![false; self.bins.len()];
    }

    fn contains(&self, element: &usize) -> bool {
        let bin_count = self.bins.len();
        for seed in self.seeds.iter() {
            if !self.bins[hash_element(element, *seed) % bin_count] {
                return false
            }
        }
        true
    }
}

#[derive(Debug, Default)]
struct Confusion {
    entries: [usize; 4],
    p2_hits: usize,
    p2_hits_target: usize,
}

impl Confusion {
    fn new() -> Self {
        Confusion {
            entries: [0; 4],  // [tn, fn, fp, tp]
            p2_hits: 0,
            p2_hits_target: 0,
        }
    }

    fn merge(self, other: Confusion) -> Confusion {
        Confusion {
            entries: [self.entries[0] + other.entries[0], self.entries[1] + other.entries[1], self.entries[2] + other.entries[2], self.entries[3] + other.entries[3]],
            p2_hits: self.p2_hits + other.p2_hits,
            p2_hits_target: self.p2_hits_target + other.p2_hits_target,
        }
    }
}

impl Display for Confusion {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "[{}, {}, {}, {}, {}] = [tp, fn, tn, fp, p2_hits_target]", self.entries[3], self.entries[1], self.entries[0], self.entries[2], self.p2_hits_target)
    }
}

fn main() {
    // Parameters
    let opt = Opt::from_args();
    let u_size: usize = 1 << opt.universe_size;
    let k: usize = 1 << opt.max_set_size;
    let u2_size: usize = k * opt.universe_client_factor;
    assert!(u2_size < u_size);
    let x2_size: usize = 1 << opt.client_set_size;
    assert!(x2_size < u2_size);
    assert!(x2_size <= k);
    let e_pow: isize = -opt.false_positives;
    let verbose: u8 = opt.verbose;
    let trials = opt.trials;
    let covert = opt.covert;

    // Setup
    let bf_params = BloomFilterParams::new(k, e_pow);
    if verbose >= 1 { println!("% Bloom filter parameters: {:?}", bf_params); }

    let u2 = 0..u2_size;
    let u1 = u2_size..u_size;

    let num_threads = num_cpus::get() - 2;
    if verbose >= 1 { println!("% Number of threads: {:?}", num_threads); }

    // Attack
    let start = Instant::now();
    let (attack_set, target) = snipe(&bf_params, u1, k, u2, x2_size, num_threads, covert, verbose);
    if verbose >= 1 { println!("% Finding attack set took: {:?}", start.elapsed()); }

    // Evaluate
    let start = Instant::now();
    let score: Confusion = (0..trials)
        .into_par_iter()
        .fold(|| Confusion::new(), |acc, _| {
            let acc2 = agg_confusion(u2_size, x2_size, &bf_params, &attack_set, target, true, acc);
            agg_confusion(u2_size, x2_size, &bf_params, &attack_set, target, false, acc2)
        })
        .reduce(|| Confusion::new(), |a, b| a.merge(b));
    if verbose >= 1 { println!("% Running {:?} trials took: {:?}", trials, start.elapsed()); }

    let accuracy = ((score.entries[3] + score.entries[0]) as f64) / (trials as f64);
    let advantage = accuracy - 1.0;
    println!("{:.5} \\\\ % {}", advantage, score)
}

fn agg_confusion(u2_size: usize, x2_size: usize, bf_params: &BloomFilterParams, attack_set: &Vec<usize>, target: usize, with_target: bool, mut acc: Confusion) -> Confusion {
    let mut rng = thread_rng();

    let mut victim_set = sample(&mut rng, u2_size, x2_size).into_vec();
    while victim_set.contains(&target) == with_target {
        victim_set = sample(&mut rng, u2_size, x2_size).into_vec();
    }

    // set the bloom filter and check if victim has the target
    let mut bloom_filter = BloomFilter::new(&bf_params);
    let mut victim_has_target = false;
    for element in &victim_set {
        bloom_filter.insert(*element);
        if *element == target {
            victim_has_target = true
        }
    }

    let mut intersection_exists: bool = false;
    for element in attack_set.iter() {
        if bloom_filter.contains(element) {
            intersection_exists = true;
            break;
        }
    }

    // If there is an intersection, we assume that the target is in the victim's set
    acc.entries[((intersection_exists as usize) << 1) + victim_has_target as usize] += 1;

    // Now let's try to check if P2 has any matches
    // We clear the bloom_filter, creating a new bloom_filter will reset the hash functions
    bloom_filter.clear();

    for &element in attack_set {
        bloom_filter.insert(element);
    }

    acc.p2_hits += victim_set.iter()
        .filter(|element| bloom_filter.contains(*element))
        .count();

    if bloom_filter.contains(&target) {
        acc.p2_hits_target += 1;
    }

    acc
}
