import streamlit as st
import pandas as pd

masses = {
    "A": 71.03711, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "K": 128.09496,
    "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841
}

mass_H2O = 18
mass_proton = 1


def digest_trypsin(sequence):
    peptides = []
    start = 0
    for i in range(len(sequence)-1):
        if sequence[i] in ("K", "R") and sequence[i+1] != "P":
            peptides.append(sequence[start:i+1])
            start = i+1
    peptides.append(sequence[start:])
    return peptides


def digest_aspN(sequence):
    peptides = []
    start = 0
    for i in range(len(sequence)):
        if sequence[i] == "D" and i != 0:  
            peptides.append(sequence[start:i])
            start = i
    peptides.append(sequence[start:])
    return peptides


def mass_calculator(peptide):
    MW = 0
    for aa in peptide:
        MW += masses[aa]
    MW += mass_H2O + mass_proton
    return MW


def y_ion_masses(peptide):
    mass_H2O = 18.010565   
    mass_proton = 1.007276 
    series = []
    current_mass = mass_H2O + mass_proton  
   
    for aa in reversed(peptide):
        current_mass += masses[aa]
        series.append(round(current_mass, 5))
    return series


def y_ion_identities(peptide):
    series = []
    fragment = ""
    for aa in reversed(peptide):
        fragment = aa + fragment
        series.append(fragment)
    return series


def translate_yions(y_ions, tolerance=0.1):
    masses_for_translation = {
    "K": 128.09496, "R": 156.10111, "N": 114.04293, "D": 115.02694,
    "C": 103.00919, "E": 129.04259, "Q": 128.05858, "G": 57.02146,
    "H": 137.05891, "I": 113.08406, "L": 113.08406, "A": 71.03711,
    "M": 131.04049, "F": 147.06841, "P": 97.05276, "S": 87.03203,
    "T": 101.04768, "W": 186.07931, "Y": 163.06333, "V": 99.06841
    }
    #I use a different AA order for this tool bcs K and Q can get confused (only 0.04 difference).
    #Fortuntely, there shouldn't be a K in the middle of a tryptic peptide so I can cheat my way
    #out of this problem
    
    y1 = round(y_ions[0], 1)
    if abs(y1 - 147.1) < 0.5:
        sequence = ["K"]
    elif abs(y1 - 175.1) < 0.5:
        sequence = ["R"]
    
    for i in range(1, len(y_ions)):
        diff = y_ions[i] - y_ions[i-1]
        
        for aa in masses_for_translation:
            mass = masses_for_translation[aa]
            if abs(diff - mass) < tolerance:
                residue = aa   
        sequence.append(residue)

    return "".join(sequence[::-1])


def findPrecursor(target, tolerance=1):

    hits = []
    aas = list(masses.keys())

    #1 aa
    for aa in aas:
        mass = masses[aa]
        if abs(target - mass) < tolerance:
            error = abs(target - mass)
            hits.append({aa : error})


    #2 aas
    for i in range(len(aas)):
        for j in range(i, len(aas)):

            aa1 = aas[i]
            aa2 = aas[j]

            mass = masses[aa1] + masses[aa2]

            if abs(target - mass) < tolerance:
                error = abs(target - mass)
                hits.append({aa1 + " + " + aa2: error})


    #3 aas
    for i in range(len(aas)):
        for j in range(i, len(aas)):
            for k in range(j, len(aas)):

                aa1 = aas[i]
                aa2 = aas[j]
                aa3 = aas[k]

                mass = (
                    masses[aa1]
                    + masses[aa2]
                    + masses[aa3]
                )

                if abs(target - mass) < tolerance:
                    error = abs(target - mass)
                    hits.append({aa1 + " + " + aa2 + " + " + aa3: error})


    return hits

st.title("Cool protein stuff")

st.subheader('Use this first box/area for proteolytic digests and y-ion series building (pretty much the things that require PROTEIN sequences)')

sequence = st.text_area('Enter a protein sequence')
protein_tools = {
    'Trypsin Digest': digest_trypsin,
    'AspN Digest': digest_aspN,
    'Y-ion Series Builder (Gives a y-ion series + masses based on a peptide)': y_ion_masses
}


protein_tool = st.selectbox('Pick a Tool (All tools will also provide masses)', protein_tools.keys())

if st.button('Digest / Analyse'):

    sequence = sequence.replace(' ', '').replace('\n', '').upper()
    if not sequence:
        st.error("Sequence cannot be empty. If you believe that the tool is being used correctly, please contact me.")
    else:
        
        validSeq = True
        for aa in sequence:
            if aa not in set("ACDEFGHIKLMNPQRSTVWY"):
                st.error("Either I programmed this incorrectly, or you have provided a sequence with letters that are NOT Amino Acids. Spaces or new-lines are okay, but try not to have your letters delimited by any other character.")
                validSeq = False
                break

        if validSeq:
            result = protein_tools[protein_tool](sequence)


            if protein_tool in ['Trypsin Digest', 'AspN Digest']:
                peptide_masses = []

                for peptide in result:
                    peptide_masses.append(mass_calculator(peptide))

                numLines = len(peptide_masses)

                st.table({
                    'Peptide': result,
                    'Mass': peptide_masses
                })

                st.download_button('Download peptide names (.txt)', "\n".join(result), f'{protein_tool.replace(' ', '_').lower()}_digested_peptides.txt', on_click='ignore')
                st.download_button('Download peptide masses (.txt)', "\n".join([str(mass) for mass in peptide_masses]), f'{protein_tool.replace(' ', '_').lower()}_peptide_masses.txt', on_click='ignore')

                df = pd.DataFrame({
                    'Peptide': result,
                    'Mass': peptide_masses
                })

                st.download_button('Download table (.csv)', df.to_csv(index=False), f'{protein_tool.replace(' ', '_').lower()}_table.csv', on_click='ignore')

                st.text(
                f"You can just download the txts and paste them into a table with {numLines} "
                "rows (not incl. header) if you dont want the ugly formatting. (As in like the black highlighting from System dark mode or csv weirdness and stuff).\n\n"
                "You could alternitavely just paste the table from here if it works out for you (I just personally prefer my own table formatting from scratch)."
            )



            else:
                if protein_tool == 'Y-ion Series Builder (Gives a y-ion series + masses based on a peptide)':


                    #I know this is a stupid control flow but I had a big idea of integrating other tools when I started
                    #and now the fastest way for me to get these done quick/dirty with accurate filenames WHILE somewhat allowing more tools later
                    #is this so deal with it. pls dont judge :pray:

                    #Note <result> is already defined as the output of the chosen tool (y ion generation in this case which is a list of y-ion masses)
                    y_ions_identity = y_ion_identities(sequence)

                    st.table({
                        'Name': [f'y{i+1}' for i in range(len(result))],
                        'Y_ion': y_ions_identity,
                        'm/z': [f"{mass:.5f}".rstrip('0').rstrip('.') for mass in result]
                        #^That's just rerouted to be more complicated bcs stupid streamlit was auto-rounding my precise calculations
                    })

                    st.download_button('Download y-ion Identities (.txt)', '\n'.join(y_ions_identity), 'y_ion_identities.txt', on_click='ignore')
                    #^I'm allowed to hardcode filename bcs of the silly control flow around tool selection (I promise the real stuff I make is cleaner and more efficient).

                    st.download_button('Download y-ion Masses (.txt)', '\n'.join([str(mass) for mass in result]), 'y_ion_masses.txt', on_click='ignore')

                    df = pd.DataFrame({
                        'name': [f'y{i+1}' for i in range(len(result))],
                        'y_ion': y_ions_identity,
                        'mass': result
                    })

                    st.download_button('Download table (.csv)', df.to_csv(index=False), 'y_ion_translation.csv', on_click='ignore')

                    st.text(
                            f"You can just download the txts and paste them into a table with {len(result)} "
                            "rows (not incl. header) if you dont want the ugly formatting. (As in like the black highlighting from System dark mode or csv weirdness and stuff).\n\n"
                            "You could alternitavely just paste the table from here if it works out for you (I just personally prefer my own table formatting from scratch)."
                        )


st.divider()

st.subheader('Use this part instead, IF working on building a protein sequence based on some given PEAKLIST\n',
             'You can leave the above part blank')


peaklist = st.text_area('Enter a peaklist here. Please try to make sure your numbers are separated by commas, and have at least one d.p of precision. Tabular or space-separated values should be accepted but comma-separated is preferred')
st.text('Example input: 147.1, 250.2, 321.3')
peaklist_tools = {'Translate y-ion Spectra to Protein Seqence': translate_yions}
st.text('Do NOT include the precursor ion mass here (full mass of the peptide). There is a tool to figure out that mass below (sometimes there is more than one AA between the last y-ion and precursor)')

peaklist_tool = st.selectbox('Pick a Tool (Theres just one here for now)', peaklist_tools.keys())
st.text('One day I might feel compelled to make a proper hub of tools, so the selection menu is just to provide some infrastructure for a potential future project')

tolerance = st.slider('Mass Tolerance (Da)', min_value=0.03, max_value=1.0, value=1.0, step=0.01)

st.markdown("""
**You should ideally keep this at 1.0**.
- From experimenting back when I made this (a year ago), 1.0 was the sweet spot that gave me a correct ion series.
- Lower values would theoretically be more accurate, but general imprecision and things beyond my comprehension introduce small errors.
- My advice would be to start at 1.0 and test the output against your own calculations to ensure correctness AND to not interpret this output as entirely robust or correct.
- I only really made this tool to see if it was possible last year, and it WAS given 1.0 tolerance for all tests, but I still don't entirely trust it, so **please use it carefully**.
- Its very janky. For instance, it couldn't accurately tell the difference between K(128.09496) and Q(128.05858), so I cheated by knowing that K won't appear in the middle of a tryptic peptide (ideally).
- Remember Isoleucine and Leucine cannot be distinguished between until you run the BLASTs yourself.
""") #Apparently st.text() doesn't like multiple lines.

if st.button('Analyse'):

    try:

        peaklist = [float(i) for i in peaklist.replace(',', " ").split()]

        if not peaklist:
            st.error("Please provide some valid peaklist in the box above. If you do believe that it is being used correctly, please contact me.")

        else:

            ion_series = peaklist_tools[peaklist_tool](peaklist, tolerance)

            st.table({'Peptide': ion_series})

    except ValueError:
        st.error("Please provide a valid peaklist containing only numbers separated by spaces or commas. If you believe that you are using the input correctly, please contact me.")

st.divider()

st.subheader('This part figures out the amino acid sequences between your last recorded y-ion mass and the precursor ion')
st.text("This should just be one amino acid (very easy to calculate), but there has been one such case in the past where multiple AAs could've fit the gap.\n")

st.text("Calculate the difference between your final recorded y-ion mass and precursor mass")
precursor_difference = st.text_area("Enter that difference here:")

if st.button('Find missing AAs'):

    if not precursor_difference.strip():
        st.error("Please enter a precursor mass difference.")

    else:

        precursor_difference_isValid = True
        try:
            precursor_difference = float(precursor_difference)
            if precursor_difference <= 0:
                raise ValueError

        except (ValueError, TypeError):
            precursor_difference_isValid = False
            st.error("Please enter a valid number, e.g. 128.09")

        if precursor_difference_isValid:
            hits = findPrecursor(precursor_difference)

            if not hits:
                st.warning("No amino acid combinations were found within the selected tolerance (1 Da). Let me know if this is raised. I don't believe it will be as of now (given appropriate input) and have not built around it yet.")

            else:
                st.table({
                    "Hit": [key for dct in hits for key in dct],
                    "Error": [val for dct in hits for val in dct.values()]
                })

            st.markdown("""
            All the hits above are possible Amino Acids that can fit into that gap.

            The error is just the difference between the mass of those chosen amino acids and the gap between precursor to final y-ion.

            **I only considered gaps of up to 3 amino acids**. Technically, four amino acids are only possible if the mass difference is:
            ```
            57.1 * 4 = ~228.4 or higher.
            ```
            ^That is the mass of glycine (smallest mass AA) x 4. Very unlikely to see that, however you should consider 4+ AAs in the gap if your precursor difference points to that outcome. If this is problematic, let me know and I can perhaps try to build a better model.

            I'd recommend running BLAST searches with those possible amino acids added to your translated spectrum to see which one fits.

            Remember hits with two amino acids (i.e A+G) can appear in either order (i.e A+G or G+A) because we only have the amino acid identities to work with (not order).
            """)


st.divider()

st.subheader("Peptide Mass Calculator")
st.text("This is just a small calculator that gives you the mass of a given peptide. Might help validate stuff when doing many calculations")
sequence_for_mass = st.text_area("Enter Sequence Here:")

if st.button("Calculate Mass"):
    sequence_for_mass = sequence_for_mass.replace(' ', '').upper()

    validSeq = True
    for aa in sequence_for_mass:
        if aa not in set("ACDEFGHIKLMNPQRSTVWY"):
            st.error("Either I programmed this incorrectly, or you have provided a sequence with letters that are NOT Amino Acids. Spaces or new-lines are okay, but try not to have your letters delimited by any other character.")
            validSeq = False
            break

    if validSeq:
        mass = mass_calculator(sequence_for_mass)

        st.table({
            'Mass (Da)': mass
        })

    

with st.sidebar:

    st.sidebar.text('NOTE: These are the masses used for all calculations')

    st.table({
            'Amino Acid': list(masses.keys()),
            'Monoisotopic Mass (Da)': list(masses.values())
        })

    st.markdown("""
    You can view the source code here:
    [Github Repository](https://github.com/Shiv276/Proteomic-tools-for-coursework)
    The code is generally sloppy and kind of inefficient, but functions well enough to provide what I needed.
    """)






