# File name: streamlit_app.py (save this as is for easy deployment)

import streamlit as st
from sympy import symbols, Eq, solve, Rational, lcm
import re
from math import gcd
from functools import reduce

# Custom CSS for a clean look
st.markdown("""
<style>
    .main {background-color: #f0f2f6;}
    .stTextInput > label {font-weight: bold;}
    .stButton > button {background-color: #4CAF50; color: white; border-radius: 8px;}
    .result {font-size: 1.2em; font-family: monospace; background-color: #e8f5e8; padding: 10px; border-radius: 5px;}
</style>
""", unsafe_allow_html=True)

def parse_formula(formula):
    atoms = {}
    charge = 0
    i = 0
    n = len(formula)
    while i < n:
        if formula[i].isupper():
            atom = formula[i]
            i += 1
            if i < n and formula[i].islower():
                atom += formula[i]
                i += 1
            num = 0
            while i < n and formula[i].isdigit():
                num = num * 10 + int(formula[i])
                i += 1
            atoms[atom] = num if num > 0 else 1
        else:
            i += 1
    
    if i < n:
        if formula[i] == '+':
            sign = 1
            i += 1
        elif formula[i] == '-':
            sign = -1
            i += 1
        else:
            sign = 0
        num_charge = 0
        while i < n and formula[i].isdigit():
            num_charge = num_charge * 10 + int(formula[i])
            i += 1
        if num_charge > 0:
            charge = sign * num_charge
        elif sign != 0:
            charge = sign
    
    return atoms, charge

def get_atoms(compounds):
    all_atoms = set()
    for c in compounds:
        atoms, _ = parse_formula(c)
        all_atoms.update(atoms.keys())
    return sorted(list(all_atoms))

def atom_count_vector(compound, atoms):
    atoms_count, _ = parse_formula(compound)
    return [atoms_count.get(a, 0) for a in atoms]

def charge_vector(compound):
    _, charge = parse_formula(compound)
    return charge

def balance_equation(reactants_str, products_str):
    reactants = [c.strip() for c in reactants_str.split('+') if c.strip()]
    products = [c.strip() for c in products_str.split('+') if c.strip()]
    
    if not reactants or not products:
        return "Invalid input: Need both reactants and products"
    
    all_compounds = reactants + products
    atoms = get_atoms(all_compounds)
    
    if not atoms:
        return "No atoms found"
    
    n_react = len(reactants)
    n_prod = len(products)
    n_total = n_react + n_prod
    
    coeffs = symbols(f'c0:{n_total}')
    c0 = coeffs[0]
    
    equations = []
    for i, atom in enumerate(atoms):
        left = sum(coeffs[j] * atom_count_vector(reactants[j], atoms)[i] for j in range(n_react))
        right = sum(coeffs[n_react + j] * atom_count_vector(products[j], atoms)[i] for j in range(n_prod))
        equations.append(Eq(left, right))
    
    left_charge = sum(coeffs[j] * charge_vector(reactants[j]) for j in range(n_react))
    right_charge = sum(coeffs[n_react + j] * charge_vector(products[j]) for j in range(n_prod))
    equations.append(Eq(left_charge, right_charge))
    
    try:
        sol_dict = solve(equations, coeffs[1:], dict=True)
        if not sol_dict:
            return "Cannot balance: No solution"
        sol = sol_dict[0]
        for c in coeffs[1:]:
            if c in sol:
                sol[c] = sol[c].subs(c0, 1)
        sol[c0] = 1
    except Exception as e:
        return f"Error solving: {e}"
    
    values = [sol[c] for c in coeffs]
    
    nums = []
    dens = []
    for v in values:
        if isinstance(v, Rational):
            nums.append(v.p)
            dens.append(v.q)
        else:
            nums.append(int(v))
            dens.append(1)
    
    lcd = reduce(lcm, dens, 1)
    
    values_int = [int(num * (lcd // den)) for num, den in zip(nums, dens)]
    
    non_zero_abs = [abs(v) for v in values_int if v != 0]
    if non_zero_abs:
        g = reduce(gcd, non_zero_abs)
        if g > 1:
            values_int = [v // g for v in values_int]
    
    first_react_idx = next((j for j in range(n_react) if values_int[j] != 0), None)
    if first_react_idx is not None:
        if values_int[first_react_idx] < 0:
            values_int = [-v for v in values_int]
    
    def format_compound(comp, coef):
        if coef == 0:
            return ""
        coef_str = f"{abs(coef)}" if abs(coef) != 1 else ""
        return f"{coef_str}{comp}" if coef_str else comp
    
    balanced_react = [format_compound(reactants[j], values_int[j]) for j in range(n_react) if values_int[j] != 0]
    balanced_prod = [format_compound(products[j], values_int[n_react + j]) for j in range(n_prod) if values_int[n_react + j] != 0]
    
    if not balanced_react or not balanced_prod:
        return "Invalid balance: No compounds"
    
    return ' + '.join(balanced_react) + ' → ' + ' + '.join(balanced_prod)

# Streamlit App
st.set_page_config(page_title="Chemical Equation Balancer", page_icon="⚗️", layout="wide")

st.title("⚗️ Chemical Equation Balancer")
st.markdown("Enter reactants and products (use + to separate compounds, include charges like Na+ or SO4 2-). The app balances atoms and charges automatically.")

col1, col2 = st.columns(2)

with col1:
    st.subheader("Reactants")
    reactants_input = st.text_input(
        "e.g., H2 + O2 or Na+ + Cl-",
        value="H2 + O2",
        help="Separate compounds with +"
    )

with col2:
    st.subheader("Products")
    products_input = st.text_input(
        "e.g., H2O or NaCl",
        value="H2O",
        help="Separate compounds with +"
    )

if st.button("Balance Equation", type="primary"):
    with st.spinner("Balancing..."):
        result = balance_equation(reactants_input, products_input)
        if "Invalid" in result or "Error" in result or "Cannot" in result:
            st.error(result)
        else:
            st.success("Balanced!")
            st.markdown(f"**{result}**", unsafe_allow_html=True)

st.markdown("---")
st.caption("Powered by SymPy • Handles redox and ionic equations via atom & charge conservation")

# Examples sidebar
with st.sidebar:
    st.header("Examples")
    st.code("H2 + O2 → H2O", language="text")
    st.code("Fe + O2 → Fe2O3", language="text")
    st.code("Na+ + Cl- → NaCl", language="text")
    st.code("KMnO4 + HCl → KCl + MnCl2 + Cl2 + H2O", language="text")
    st.markdown("[Deploy Guide](#)")
