# File name: streamlit_app.py (save this as is for easy deployment)

import streamlit as st
from sympy import symbols, Eq, solve, Rational, lcm
import re
from math import gcd
from functools import reduce
from collections import OrderedDict

# Custom CSS for a clean look
st.markdown("""
<style>
    .main {background-color: #f0f2f6;}
    .stTextInput > label {font-weight: bold;}
    .stButton > button {background-color: #4CAF50; color: white; border-radius: 8px;}
    .result {font-size: 1.2em; font-family: monospace; background-color: #e8f5e8; padding: 10px; border-radius: 5px;}
    .ox-table {background-color: #f9f9f9; padding: 10px; border-radius: 5px;}
    .suggestion {background-color: #fff3cd; padding: 10px; border-radius: 5px; margin: 5px 0;}
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
    
    # Handle charge at the end
    if i < n and (formula[i] in '+-'):
        sign = 1 if formula[i] == '+' else -1
        i += 1
        num_charge = 0
        while i < n and formula[i].isdigit():
            num_charge = num_charge * 10 + int(formula[i])
            i += 1
        charge = sign * (num_charge if num_charge > 0 else 1)
    
    return atoms, charge

def calculate_oxidation_states(formula):
    atoms, charge = parse_formula(formula.strip())
    if not atoms:
        return "Invalid formula"
    
    # Default oxidation states (simple rules)
    ox_defaults = {
        'H': 1, 'O': -2, 'F': -1, 'Cl': -1, 'Br': -1, 'I': -1,
        'Li': 1, 'Na': 1, 'K': 1, 'Rb': 1, 'Cs': 1, 'Fr': 1,
        'Be': 2, 'Mg': 2, 'Ca': 2, 'Sr': 2, 'Ba': 2, 'Ra': 2,
        # Add more as needed (e.g., Al:3, etc.)
    }
    
    known_sum = 0
    known_elements = []
    unknown_elements = []
    
    for elem in atoms:
        if elem in ox_defaults:
            ox = ox_defaults[elem]
            known_sum += ox * atoms[elem]
            known_elements.append(elem)
        else:
            unknown_elements.append(elem)
    
    ox_states = {}
    
    if len(unknown_elements) > 1:
        return "Too many unknown elements; specify more rules or use a simpler formula."
    
    if unknown_elements:
        unknown = unknown_elements[0]
        count_u = atoms[unknown]
        if count_u == 0:
            return "No atoms found for unknown element."
        ox_u = (charge - known_sum) / count_u
        if ox_u != int(ox_u):
            return f"Non-integer oxidation state for {unknown}: {ox_u}. May need special rules (e.g., peroxides)."
        ox_states[unknown] = int(ox_u)
    else:
        # All known
        if known_sum != charge:
            return f"Inconsistent: Sum of known ox states ({known_sum}) != charge ({charge}). Special case?"
        for elem in known_elements:
            ox_states[elem] = ox_defaults[elem]
    
    # Add knowns
    for elem in known_elements:
        ox_states[elem] = ox_defaults[elem]
    
    # Sort by appearance order (simple alphabetical for now)
    return dict(sorted(ox_states.items()))

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

def suggest_missing_compounds(reactants_str, products_str):
    reactants = [c.strip() for c in reactants_str.split('+') if c.strip()]
    products = [c.strip() for c in products_str.split('+') if c.strip()]
    
    if not reactants or not products:
        return ["Invalid input: Need both sides."]
    
    all_atoms = get_atoms(reactants + products)
    
    # Assume coef=1 for all, compute total counts
    left_count = {a: 0 for a in all_atoms}
    for r in reactants:
        atoms, _ = parse_formula(r)
        for a, c in atoms.items():
            left_count[a] += c
    
    right_count = {a: 0 for a in all_atoms}
    for p in products:
        atoms, _ = parse_formula(p)
        for a, c in atoms.items():
            right_count[a] += c
    
    delta = {a: left_count[a] - right_count[a] for a in all_atoms}
    
    left_charge_total = sum(charge_vector(r) for r in reactants)
    right_charge_total = sum(charge_vector(p) for p in products)
    charge_delta = left_charge_total - right_charge_total
    
    suggestions = []
    
    # Handle charge imbalance first
    if charge_delta != 0:
        if charge_delta > 0:
            suggestions.append("Charge excess on left (+). Suggestion: Add anions (e.g., Cl⁻ or OH⁻) to products.")
        else:
            suggestions.append("Charge excess on right (-). Suggestion: Add cations (e.g., H⁺ or Na⁺) to products.")
    
    # Atom imbalances
    for atom in sorted(delta, key=lambda x: abs(delta[x]), reverse=True):  # Largest first
        d = delta[atom]
        if d == 0:
            continue
        
        if abs(d) < 1:  # Ignore tiny
            continue
        
        if atom == 'H' and d > 0:
            h_coef = d // 2 + (1 if d % 2 else 0)  # Round up for even
            suggestions.append(f"H excess on left ({d}). Suggestion: Add {h_coef} H₂O to products (common in aqueous reactions).")
        elif atom == 'H' and d < 0:
            h_deficit_left = abs(d)
            suggestions.append(f"H deficit on left ({h_deficit_left}). Suggestion: Add {h_deficit_left // 2 + (1 if h_deficit_left % 2 else 0)} H₂ to reactants.")
        
        elif atom == 'O' and d > 0:
            o_coef = d // 2 + (1 if d % 2 else 0)
            suggestions.append(f"O excess on left ({d}). Suggestion: Add {o_coef} O₂ to products (e.g., combustion).")
        elif atom == 'O' and d < 0:
            suggestions.append(f"O deficit on left ({abs(d)}). Suggestion: Add {abs(d)} O from H₂O or CO₂ to reactants.")
        
        elif atom == 'C' and d > 0:
            c_coef = d
            suggestions.append(f"C excess on left ({d}). Suggestion: Add {c_coef} CO₂ to products (common in combustion).")
        elif atom == 'C' and d < 0:
            suggestions.append(f"C deficit on left ({abs(d)}). Suggestion: Add {abs(d)} C from CH₄ or similar to reactants.")
        
        # Add more common atoms if needed, e.g., N -> N2, etc.
        elif atom == 'N' and d > 0:
            n_coef = d // 2
            suggestions.append(f"N excess on left ({d}). Suggestion: Add {n_coef} N₂ to products.")
    
    if not suggestions:
        suggestions.append("No obvious imbalances detected. Check formulas or reaction type.")
    
    return suggestions[:3]  # Limit to top 3 for brevity

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
            suggestions = suggest_missing_compounds(reactants_str, products_str)
            return f"Cannot balance: No solution.\n\n**Suggestions for missing compounds:**\n" + '\n'.join([f"• {s}" for s in suggestions])
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
        coef_str = f"{abs(coef)} " if abs(coef) != 1 else ""
        return f"{coef_str}{comp}" if coef_str else f"{comp}"
    
    balanced_react = [format_compound(reactants[j], values_int[j]) for j in range(n_react) if values_int[j] != 0]
    balanced_prod = [format_compound(products[j], values_int[n_react + j]) for j in range(n_prod) if values_int[n_react + j] != 0]
    
    if not balanced_react or not balanced_prod:
        return "Invalid balance: No compounds"
    
    return ' + '.join(balanced_react) + ' → ' + ' + '.join(balanced_prod)

# Streamlit App
st.set_page_config(page_title="Chemical Tools: Balancer & Oxidation States", page_icon="⚗️", layout="wide")

st.title("⚗️ Chemical Equation Tools")
st.markdown("Balance equations or calculate oxidation states. Handles atoms, charges, and simple rules. Now suggests missing compounds like H₂O if unbalanced!")

tab1, tab2 = st.tabs(["🧪 Balance Equation", "🔢 Oxidation States"])

with tab1:
    st.markdown("Enter reactants and products (use + to separate compounds, include charges like Na+ or SO4 2-). If unbalanced, get suggestions for missing parts!")
    
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
            if "Cannot balance" in result or "Error" in result:
                st.error("**Balancing failed.**")
                st.markdown(result.replace('\n', '<br>'), unsafe_allow_html=True)
            else:
                st.success("Balanced!")
                st.markdown(f"**{result}**", unsafe_allow_html=True)

with tab2:
    st.markdown("Enter a chemical formula to calculate oxidation states using standard rules (H=+1, O=-2, etc.). Handles charges like SO4 2-. Note: Simple rules; doesn't handle peroxides or complex cases.")
    
    formula_input = st.text_input(
        "Formula",
        value="KMnO4",
        help="e.g., H2O, CH4, KMnO4, SO4 2-"
    )
    
    if st.button("Calculate Oxidation States", type="primary"):
        with st.spinner("Calculating..."):
            result = calculate_oxidation_states(formula_input)
            if isinstance(result, str):
                if "Invalid" in result or "Too many" in result or "Inconsistent" in result:
                    st.error(result)
                else:
                    st.warning(result)
            else:
                st.success("Oxidation states calculated!")
                # Display as table
                df_data = [{"Element": elem, "Oxidation State": ox} for elem, ox in result.items()]
                st.dataframe(df_data, use_container_width=True, hide_index=True)
                # Also show sum check
                atoms, charge = parse_formula(formula_input)
                total = sum(ox * atoms.get(elem, 0) for elem, ox in result.items())
                st.info(f"Total sum: {total} (should equal charge: {charge})")

# Examples sidebar
with st.sidebar:
    st.header("Examples")
    st.markdown("**Balancer:**")
    st.code("H2 + O2 → H2O", language="text")
    st.code("CH4 + O2 → CO2 (try this - will suggest missing H2O!)", language="text")
    st.code("Na+ + Cl- → NaCl", language="text")
    st.markdown("**Ox States:**")
    st.code("KMnO4 → K:+1, Mn:+7, O:-2", language="text")
    st.code("CH4 → C:-4, H:+1", language="text")
    st.markdown("---")
    st.caption("Powered by SymPy • Simple ox rules")

# Footer
st.markdown("---")
st.caption("Deployed with Streamlit • Suggestions use heuristics for common cases like H₂O in aqueous/combustion")
