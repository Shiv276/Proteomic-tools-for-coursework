aa_masses = {
    "K": 128.09496, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "A": 71.03711,
    "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841
}

def translate_yions(y_ions):
    y1 = round(y_ions[0], 1)
    if abs(y1 - 147.1) < 0.5:
        sequence = ["K"]
    elif abs(y1 - 175.1) < 0.5:
        sequence = ["R"]
    
    for i in range(1, len(y_ions)):
        diff = y_ions[i] - y_ions[i-1]
        
        for aa in aa_masses:
            mass = aa_masses[aa]
            if abs(diff - mass) < 0.1:
                residue = aa   
        sequence.append(residue)

    return "".join(sequence[::-1])

y_ions_1_test = [147.1, 303.18, 417.3, 545.37, 674.46, 745.5, 
858.57, 986.65, 1115.71, 1278.76]

y_ions_2_test = [147.07, 276.15, 347.19, 446.26, 533.32, 604.37, 
703.40, 760.44, 889.51, 988.58, 1087.61, 1144.65]

y_ions_3_test = [175.09, 338.19, 494.30, 609.39, 708.41, 
765.45, 852.53, 951.59, 1066.61, 1165.65]

y_ions_1 = translate_yions(y_ions_1_test)
y_ions_2 = translate_yions(y_ions_2_test) 
y_ions_3 = translate_yions(y_ions_3_test)


with open("y_ion_series_1.txt", "w") as f:
    for fragments in y_ions_1:
        f.write(str(fragments))

with open("y_ion_series_2.txt", "w") as f:
    for fragments in y_ions_2:
        f.write(str(fragments))

with open("y_ion_series_3.txt", "w") as f:
    for fragments in y_ions_3:

        f.write(str(fragments))
