def remove_and_create_new_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Split the line into columns
            columns = line.strip().split()

            # Check if there are at least three columns
            if len(columns) >= 3:
                # Keep only the first two columns
                new_line = ' '.join(columns[:2]) + '\n'

                # Write the modified line to the new file
                outfile.write(new_line)

# Replace 'input.txt' with the name of your input file
# Replace 'output.txt' with the desired name for the new file
remove_and_create_new_file('133Cs_only_exp.txt', 'output.txt')

# Initial guess for parameters
# initial_params = [95.83, 52169.0, 1e-13 , 1.0 , 2019.0]
# initial_params = [80,52169.65,2001,22.438,5.5e-13]
#
# initial_params = [-0.8564708287207337,110608.53651787751,2006.0272437830115,19.97941955508578, 43.57813758756277]