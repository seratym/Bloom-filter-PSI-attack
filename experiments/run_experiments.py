import subprocess

executable = "../target/release/bloom-goes-boom"
trials = 1000000
universe_size = 32

options_priority = [
    (5, 8, 2),
    (5, 8, 3),
    (5, 8, 4),
    # (5, 12, 2),
    # (5, 12, 3),
    # (5, 12, 4),
    # (5, 16, 2),
    # (5, 16, 3),
    # (5, 16, 4),
    (10, 8, 2),
    (10, 8, 3),
    (10, 8, 4),
    (10, 12, 2),
    # (10, 12, 3),
    # (10, 12, 4),
    (10, 16, 2),
    # (10, 16, 3),
    # (10, 16, 4),
    (20, 8, 2),
    (20, 8, 3),
    (20, 8, 4),
    (20, 12, 2),
    # (20, 12, 3),
    # (20, 12, 4),
    (20, 16, 2),
    # (20, 16, 3),
    # (20, 16, 4),
    (30, 8, 2),
    (30, 8, 3),
    (30, 8, 4),
    (30, 12, 2),
    # (30, 12, 3),
    # (30, 12, 4),
    (30, 16, 2),
    # (30, 16, 3),
    # (30, 16, 4),
]

# Function to run the bloom-goes-boom binary with given parameters
def run_bloom(false_positive, max_set_size, universe_client_factor):
    # Construct the command as a list of arguments
    flags = [
        f"--universe-size={universe_size}",
        f"--universe-client-factor={universe_client_factor}",
        f"--client-set-size={max_set_size}",
        f"--max-set-size={max_set_size}",
        f"--false-positives={false_positive}",
        f"--trials={trials}",
        "-v",
    ]
    cmd = [executable] + flags
    print(f"% {' '.join(cmd)}")
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    print(f"{false_positive} & {max_set_size} & {universe_client_factor} &  % fp k U2_fact")
    for line in process.stdout:
        print(line, end="")  # Print each line of output as it's generated
    for line in process.stderr:
        print(line, end="")  # Print each line of output as it's generated
    process.wait()

for i, (false_positive, max_set_size, universe_client_factor) in enumerate(options_priority):
    run_bloom(false_positive, max_set_size, universe_client_factor)
