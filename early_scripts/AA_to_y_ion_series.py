aa_masses = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841
}

mass_H2O = 18.010565   
mass_proton = 1.007276 

def y_ion_series(peptide):
    series = []
    current_mass = mass_H2O + mass_proton  
   
    for aa in reversed(peptide):
        current_mass += aa_masses[aa]
        series.append(round(current_mass, 5))
    return series


peptide = "RALISSTGPYWADKR"
y_ions = y_ion_series(peptide)


with open("y_ion_series_masses.txt", "w") as f:
    for fragments in y_ions:
        f.write(str(fragments) + "\n")


def y_ion_identities(peptide):
    
    series = []
    fragment = ""
    for aa in reversed(peptide):
        fragment = aa + fragment
        series.append(fragment)
    return series

y_identities = y_ion_identities(peptide)

with open("y_ion_series_identities.txt", "w") as f:
    for frag in y_identities:
        f.write(frag + "\n")
