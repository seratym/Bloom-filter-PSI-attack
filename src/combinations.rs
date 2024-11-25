use fftw::array::AlignedVec;
use fftw::plan::*;
use fftw::types::*;
use num::integer::binomial;
use num::BigUint;
use num::FromPrimitive;
use rand::thread_rng;
use rand::Rng;
use std::cmp::min;
use std::fmt;
use std::fmt::Formatter;
use std::ops::Range;


struct FftPlans {
    n_pow2: usize,
    forward: C2CPlan64,
    backward: C2CPlan64,
    uuid: [u8; 16],
}

impl FftPlans {
    pub fn new(n_pow2: usize) -> FftPlans {
        let forward = C2CPlan::aligned(&[n_pow2], Sign::Forward, Flag::PATIENT | Flag::DESTROYINPUT).unwrap();
        let backward = C2CPlan::aligned(&[n_pow2], Sign::Backward, Flag::PATIENT | Flag::DESTROYINPUT).unwrap();

        FftPlans {
            n_pow2,
            forward,
            backward,
            uuid: thread_rng().gen(),
        }
    }
}

struct Polynomial {
    coefficients: AlignedVec<c64>,
    uuid: [u8; 16],
}
impl fmt::Debug for Polynomial {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.coefficients.to_vec())
    }
}

impl Polynomial {
    pub fn from_slice(coefficents: &[usize], fft_plans: &mut FftPlans) -> Self {
        let mut aligned_vec = AlignedVec::new(fft_plans.n_pow2);

        for i in 0..coefficents.len() {
            aligned_vec[i] = c64::new(coefficents[i] as f64, 0.0);
        }

        Self::new(aligned_vec, fft_plans)
    }

    pub fn new(mut aligned_vec: AlignedVec<c64>, fft_plans: &mut FftPlans) -> Self {
        let mut fft_coefficients = AlignedVec::new(fft_plans.n_pow2);
        fft_plans.forward.c2c(&mut aligned_vec, &mut fft_coefficients).unwrap();

        Polynomial {
            coefficients: fft_coefficients,
            uuid: fft_plans.uuid,
        }
    }

    pub fn multiply(mut self, other: Polynomial) -> Polynomial {
        assert_eq!(self.uuid, other.uuid);

        for i in 0..self.coefficients.len() {
            self.coefficients[i] *= other.coefficients[i];
        }

        Polynomial {
            coefficients: self.coefficients,
            uuid: self.uuid,
        }
    }

    pub fn get_coefficients(mut self, indices: Range<usize>, fft_plans: &mut FftPlans) -> Vec<usize> {
        assert_eq!(self.uuid, fft_plans.uuid);
        let mut aligned_vec = AlignedVec::new(fft_plans.n_pow2);

        fft_plans.backward.c2c(&mut self.coefficients, &mut aligned_vec).unwrap();

        aligned_vec[indices].iter().map(|x| (x.re / fft_plans.n_pow2 as f64).round() as usize).collect()
    }
}

fn generate_bin_polynomial(elements_in_bin: usize, fft_plans: &mut FftPlans) -> Polynomial {
    let mut aligned_vec = AlignedVec::new(fft_plans.n_pow2);

    for j in 1..=elements_in_bin {
        aligned_vec[j - 1] = c64::new(binomial(elements_in_bin, j) as f64, 0.0);
    }

    Polynomial::new(aligned_vec, fft_plans)
}

fn product(polynomials: Vec<Polynomial>) -> Polynomial {
    polynomials.into_iter().reduce(|a, b| a.multiply(b)).unwrap()
}

fn compute_appropriate_pow2(elements_in_bins: &[usize]) -> usize {
    let sum: usize = elements_in_bins.iter().sum();
    (sum - elements_in_bins.len() + 1 + 1).next_power_of_two()
}

pub fn compute_unique(elements_in_bins: &[usize], page_size: usize, universe_size: usize) -> BigUint {
    if elements_in_bins.is_empty() {
        let r = elements_in_bins.len();
        return binomial(BigUint::from_usize(universe_size).unwrap(), BigUint::from_usize(page_size - r).unwrap())
    }

    // LHS
    let n = compute_appropriate_pow2(elements_in_bins);
    let fft_plans = &mut FftPlans::new(n);
    let polynomials = elements_in_bins.iter().map(|x| generate_bin_polynomial(*x, fft_plans)).collect();
    let out = product(polynomials);

    // RHS
    let r = elements_in_bins.len();
    let sum: usize = elements_in_bins.iter().sum();
    let s = sum - r + 1;
    let smallest = min(s, page_size - r);
    let universe_prime_size = universe_size - sum;
    let lhs = out.get_coefficients(0..(smallest + 1), fft_plans);
    lhs.into_iter().enumerate().map(|(j, a)| a * binomial(BigUint::from_usize(universe_prime_size).unwrap(), BigUint::from_usize(page_size - r - j).unwrap())).sum()
}



#[test]
fn test_easy_poly() {
    let easy_options: Vec<(usize, usize, Vec<usize>, &[u8])> = vec![
        (2, 2, vec![1, 1], b"1"),
        (11, 4, vec![2, 2, 3, 4], b"48"),
        (10, 4, vec![2, 2, 3], b"60"),
        (8, 4, vec![2, 2], b"41"),
        (8, 5, vec![2, 2], b"44"),
        (6, 3, vec![1, 2], b"7"),
        (6, 3, vec![3, 2], b"15"),
        (11, 4, vec![2, 2, 3, 3], b"36"),
        (11, 4, vec![2, 2, 3, 4], b"48"),
        (10, 4, vec![2, 2, 3, 3], b"36"),
    ];

    for option in easy_options.iter() {
        let (universe_size, page_size, bins, expected_result) = option;
        let result = compute_unique(bins, *page_size, *universe_size);
        assert_eq!(result, BigUint::parse_bytes(*expected_result, 10).unwrap(), "Failed with: {:?}", option);
    }
}

#[test]
fn test_hard_poly() {
    let hard_options: Vec<(usize, usize, Vec<usize>, &[u8])> = vec![
        (1024, 512, vec![1, 1, 3, 1, 1, 3, 1, 2, 1, 1, 4, 1, 2, 1, 3, 1, 2, 1, 1, 1],
         b"127745056341237764056345943210246739431214534442373457326761996111054660211302977467973680840811954324644997330243445168657301827328550016587576171036184140641980322809527743907896223344207125636066303731426559082109021141981238090860102733907929416366309897810232083955982405818052607169345566690872835"),
        (1024, 512, vec![2, 2, 2, 3, 3, 3, 1, 1, 3, 2, 1, 4, 2, 2, 4, 2, 2, 4, 1, 2],
         b"9492954612413596156955847391304088506275941943039651639226791889752393767490689810902185550963694881894176392707411301899765985021889266441782180531331666392663061513637269584467155338674521230232738771373105487332293502022927409979618954415625907699572396786371650969529511977263289585989914366259133235"),
        (1024, 512, vec![1, 4, 1, 1, 2, 1, 1, 4, 2, 2, 3, 1, 2, 2, 1, 2, 1, 3, 1, 2],
         b"709482249914585558232544718858129195339929043021403250280065799013715331152825012114023804835287245332801676980240514679904635612722752644988412001780257626065648252414493394646417479237402423091978389179365593024496850838596375656313160806180101293905953038002849765639335274192017341017141390068073875"),
        (1024, 512, vec![3, 2, 1, 2, 4, 2, 3, 1, 1, 2, 1, 1, 3, 1, 3, 3, 2, 1, 1, 3],
         b"1597805479511332789652238750304741927232641054368971687199812314554208040114497163989302798306093558967741347897932525385397614813521138743128799543307589590391155570394799413784403194457224143015321563337338139303773905263742593973816236707568115589315500831054927209079344962841030895621129978647209235"),
        (1024, 512, vec![5, 3, 1, 2, 1, 3, 2, 1, 2, 3, 2, 3, 1, 1, 1, 1, 2, 1, 1, 1],
         b"530939980700543440725341932341550311898187450130038464716850555741495004696899518053254461309404369891861686293967348976909237487567491722295131070822949429184063250380072995472945201974428300879146879433572815519688003920215684450042201364443797105927075902325656603767378955442588328537318605776319875"),
        (1024, 512, vec![1, 2, 1, 4, 1, 2, 1, 2, 2, 2, 2, 4, 1, 1, 1, 2, 1, 1, 1, 3],
         b"401885752862379627267581070002816412765405200013960730083080068904746445200233957026668536823662637856586925622111952323348518074729227074792566818930563019456551594953577645457829350332090302007834874301580969695335644455756134954638778970864275682674463926848910254078018482574294118648108502825848595"),
        (1024, 512, vec![1, 2, 3, 1, 1, 1, 1, 2, 1, 1, 1, 3, 1, 1, 2, 2, 3, 1, 1, 1],
         b"101524512997979913461310879321387835996584217284253459522491897243302525226633369957548610468277905652497177842268744738605451711337169385173437065633752709119577370859135265357334202528947196124997712354674783657817638044749480311885804699185851822120951718514905885029284111680452784417222712691429955"),
    ];

    for option in hard_options.iter() {
        let (universe_size, page_size, bins, expected_result) = option;
        let result = compute_unique(bins, *page_size, *universe_size);
        assert_eq!(result, BigUint::parse_bytes(*expected_result, 10).unwrap(), "Failed with: {:?}", option);
    }
}