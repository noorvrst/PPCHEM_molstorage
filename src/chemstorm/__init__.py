# src/chemstorm/__init__.py
"""A simple tool to help students analyze 
and understand the safety risks of different molecules 
and their proper storage."""

# Define the version explicitly
__version__ = "0.1.0"

# Your current imports
from .chemstorm import *  # or your specific imports

# The launch_app function
def launch_app():
    """
    Launches the Streamlit user interface for molstorage.
    This function starts a local server and opens the application in the browser.
    """
    import streamlit.web.cli as stcli
    import sys
    import os
    
    # Find the path to app.py which is in the same folder as __init__.py
    dirname = os.path.dirname(__file__)
    app_path = os.path.join(dirname, 'app.py')
    
    # Configure the arguments for streamlit
    sys.argv = ["streamlit", "run", app_path]
    
    # Launch the application
    stcli.main()