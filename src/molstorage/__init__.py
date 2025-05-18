"""Storage of hazardous chemicals based on their pictograms, hazards statements, acid/base classification, state"""

# src/molstorage/__init__.py

from .molstorage import *  

# Add this function
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