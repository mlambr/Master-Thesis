import random

def split_file_fixed(input_file, train_file, val_file, test_file, train_count=800, val_count=200, test_count=200, seed=42):
    """
    Splits a text file into fixed numbers of training, validation, and test lines.

    Args:
        input_file (str): Path to the input file.
        train_file (str): Path to the training file output.
        val_file (str): Path to the validation file output.
        test_file (str): Path to the test file output.
        train_count (int): Number of lines to use for training.
        val_count (int): Number of lines to use for validation.
        test_count (int): Number of lines to use for testing.
        seed (int): Random seed for reproducibility (default: 42).
    """
    random.seed(seed)

    # Read all lines from the input file
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    # Ensure there are enough lines for the split
    total_lines = train_count + val_count + test_count
    if len(lines) < total_lines:
        raise ValueError(f"Not enough lines in the file. Required: {total_lines}, Found: {len(lines)}")

    # Randomly sample lines for training, validation, and testing
    selected_lines = random.sample(lines, total_lines)
    train_lines = selected_lines[:train_count]
    val_lines = selected_lines[train_count:train_count + val_count]
    test_lines = selected_lines[train_count + val_count:]

    # Write the training, validation, and test files
    with open(train_file, 'w') as train_out:
        train_out.writelines(train_lines)

    with open(val_file, 'w') as val_out:
        val_out.writelines(val_lines)

    with open(test_file, 'w') as test_out:
        test_out.writelines(test_lines)

    print(f"Split completed: {len(train_lines)} lines in training, {len(val_lines)} lines in validation, {len(test_lines)} lines in testing.")

# Paths to the files
input_file_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/Data/chromosome8.txt'
train_file_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/Data/training.txt'
val_file_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/Data/validation.txt'
test_file_path = '/cluster/work/medinfmk/cfDNA-Snyder/results/Maria/Data/testing.txt'

# Perform the split
split_file_fixed(
    input_file=input_file_path,
    train_file=train_file_path,
    val_file=val_file_path,
    test_file=test_file_path,
    train_count=800,
    val_count=200,
    test_count=400
)
