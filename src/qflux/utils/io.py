import requests
import os
from pathlib import Path

def load_dataset(ds_path: Path) -> None:
    """Loads a dataset from the specified path.

    This function is intended to load a dataset from the given file path.
    Currently, it is not implemented and always raises a NotImplementedError.

    Args:
        ds_path: The path to the dataset file.

    Raises:
        NotImplementedError: Always raised as the function is not yet implemented.
    """
    raise NotImplementedError("Not implemented")

def download_github_directory(repo_owner, repo_name, dir_path, local_dir):

    """
    Recursively download a directory from a GitHub repository
    """

    #-- GitHub API URL for the directory contents
    url = f"https://api.github.com/repos/{repo_owner}/{repo_name}/contents/{dir_path}"
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Error: {response.status_code} - {response.json().get('message', 'Unknown error')}")
        return
    
    contents = response.json()
    
    #-- Create the local directory if it doesn't exist
    os.makedirs(local_dir, exist_ok=True)
    
    for item in contents:

        if item['type'] == 'file':
            # Download the file
            file_url = item['download_url']
            file_name = item['name']
            file_path = os.path.join(local_dir, file_name)
            
            print(f"Downloading: {file_path}")
            file_response = requests.get(file_url)
            with open(file_path, 'wb') as f:
                f.write(file_response.content)
                
        elif item['type'] == 'dir':
            # Recursively download subdirectory
            subdir_path = item['path']
            subdir_local = os.path.join(local_dir, item['name'])
            print(f"===> Directory: {subdir_path}")
            download_github_directory(repo_owner, repo_name, subdir_path, subdir_local)
