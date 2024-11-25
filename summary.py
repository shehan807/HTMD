import os

def summarize_directory(root_dir, file_handle, indent=0):
    """
    Recursively summarizes the directory structure of root_dir and writes it to a file.
    
    Parameters:
        root_dir (str): The root directory to summarize.
        file_handle (file): The file handle to write the summary.
        indent (int): The level of indentation (for recursive calls).
 
    Made with ChatGPT 4o, 09/09/2024
    """
    # List all files and directories in the current root_dir
    with os.scandir(root_dir) as entries:
        for entry in entries:
            if entry.is_dir(follow_symlinks=False):
                # Write the directory name with indentation
                file_handle.write('    ' * indent + f"[DIR] {entry.name}\n")
                # Recursively summarize the subdirectory
                summarize_directory(entry.path, file_handle, indent + 1)
            elif entry.is_file(follow_symlinks=False):
                # Write the file name with indentation
                file_handle.write('    ' * indent + f"[FILE] {entry.name}\n")

if __name__ == "__main__":
    # Define the root directory (change this to the path of your main directory)
    root_directory = './'
    
    # Define the output file name
    output_file = 'summary.txt'
    
    # Check if the root directory exists
    if os.path.exists(root_directory):
        # Open the output file in write mode
        with open(output_file, 'w') as file:
            file.write(f"Directory summary for: {root_directory}\n")
            file.write('-' * 50 + '\n')
            summarize_directory(root_directory, file)
        
        print(f"Directory structure summary saved to {output_file}")
    else:
        print(f"Directory {root_directory} does not exist.")

