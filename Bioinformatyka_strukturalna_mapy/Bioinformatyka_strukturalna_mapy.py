from Bio import PDB
import matplotlib.pyplot as plt
import sys
import os


def calculate_contact_map(structure, threshold=8.0):
    contact_map = []

    # Wybierz atom CA jako reprezentanta kazdej reszty
    atoms_ca = [atom for atom in structure.get_atoms() if atom.id == 'CA']

    for atom1 in atoms_ca:
        contacts = []
        for atom2 in atoms_ca:
            if atom1 != atom2:
                # Oblicz odleglosc miedzy atomami CA
                distance = atom1 - atom2
                if distance < threshold:
                    contacts.append(1)
                else:
                    contacts.append(0)
        contact_map.append(contacts)

    return contact_map

def plot_contact_map(contact_map):
    plt.imshow(contact_map, cmap='binary', interpolation='none')
    plt.xlabel('Residues')
    plt.ylabel('Residues')
    plt.title('Contact Map')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <pdb_file>")
        sys.exit(1)

    pdb_file = sys.argv[1]
    filename = os.path.basename(pdb_file)
    pdb_id = filename[0:4]
    
    # Wczytaj strukture bialka (np. 1HHB)
    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure(pdb_id, filename)

    # Oblicz mape kontaktow
    contact_map = calculate_contact_map(structure)

    # Wygeneruj wykres
    plot_contact_map(contact_map)

