![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">ChemStorM - Chemical Storage Manager</h1>

#### A simple tool to help students analyze and understand the safety risks of different molecules and their proper storage.


## How to Install It
1. Create and activate a new environment (another name can be given).
<pre>
<code>
conda create -n chemstorm_env python=3.10 
</code>
</pre>

<pre>
<code>
conda activate chemstorm_env
</code>
</pre>


2. Pip install the package by copying this into your command line.
<pre>
<code>
(chemstorm_env) $ pip install git+https://github.com/noorvrst/chemstorm.git
</code>
</pre>

## How to Use the Streamlit Web Application
After installing, create a new file and make sure the <code>chemstorm_env</code> environment is activated. Then, run the following code to launch the Streamlit application in your browser.

```python
import chemstorm

chemstorm.launch_app()
```

Here’s an overview of the interface.

![Illustration](/assets/chemstorm_app_readme.png)

## Run Tests and Coverage
<pre>
<code>
(conda_env) $ pip install tox
(conda_env) $ tox
</code>
</pre>

## Developpers
- Verstraete Noor, student in chemistry and chemical engineering at EPFL (Switzerland), email: noor.verstraete@epfl.ch
   
[![jhc github](https://img.shields.io/badge/GitHub-noorvrst-181717.svg?style=flat&logo=github)](https://github.com/noorvrst)

- Mercier Charlotte, student in chemistry and chemical engineering at EPFL (Switzerland), email: charlotte.mercier@epfl.ch
  
[![jhc github](https://img.shields.io/badge/GitHub-chacha333%20create-181717.svg?style=flat&logo=github)](https://github.com/chacha333-create)

- Delacou Daphné, student in chemistry and chemical engineering at EPFL (Switzerland), email: daphne.delacou@epfl.ch
  
[![jhc github](https://img.shields.io/badge/GitHub-ddelacou-181717.svg?style=flat&logo=github)](https://github.com/ddelacou)
    
