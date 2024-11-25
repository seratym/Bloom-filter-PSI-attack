import itertools

# Define the ranges for the parameters
false_positives = [5, 10, 20, 30]
max_set_sizes = [8, 12, 16]
client_universe_factors = [2,3,4]

# Create all combinations of the parameter values
parameter_combinations = itertools.product(false_positives, max_set_sizes, client_universe_factors)

if __name__ == "__main__":
    for (fp, k, f) in parameter_combinations:
        print(f"({fp}, {k}, {f}),")
