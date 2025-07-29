 # Your sector data
sector_indices =  [86, 90, 91, 93, 94, 96, 101, 103, 104, 105, 106, 107, 108, 109, 110, 112, 113, 114, 117, 118, 119, 121, 125, 126, 128, 131, 132, 133, 140, 143, 144, 145, 146, 147, 149, 150, 153, 154, 155, 156, 158, 159, 162, 165, 166, 167, 169, 176, 187, 188, 190, 191, 194, 204, 208, 209, 211, 212, 213, 215, 218, 231, 232, 234, 236, 237, 238, 239, 240, 243, 244, 246, 247, 248, 249, 250, 251, 252, 254, 255, 256, 257, 258, 260, 261, 262, 263, 264, 266, 267, 271, 272, 273, 274, 275, 276, 277, 278, 279, 280, 281, 283, 284, 285, 287, 288, 289, 290, 291, 292, 293, 294, 295, 296, 298, 299, 300, 301, 304, 308, 309, 310, 311, 313, 315, 316, 317, 318, 319, 320, 321, 322, 323, 325, 327, 328, 330, 332, 334, 336, 337, 338, 339, 340, 345, 346, 348, 350, 356, 357, 358, 360, 362, 363, 364, 365, 366, 367, 368, 369, 370, 372, 373, 374, 384, 385, 388, 390, 392, 393, 394, 396, 398, 399, 401, 427, 428, 430, 462, 465, 466, 467, 471, 480, 482, 485, 486, 487, 488, 500, 502, 503, 504, 518, 520, 539, 540, 541, 563, 565, 566, 583, 584, 585, 586, 588, 589, 592, 593, 594, 595, 598, 600, 603, 610, 612, 616, 617, 618, 619, 620, 621, 622, 623, 625, 626, 628, 637, 638, 639, 640, 641, 643, 644, 647, 649, 650, 651, 652, 653, 655, 656, 657, 658, 660, 661, 662, 663, 664, 666, 669, 670, 678, 679, 680, 681, 683, 685, 686, 687, 688, 689, 690, 691, 692, 693, 694, 695, 696, 697, 699, 700, 703, 704, 708, 709, 710, 712, 713, 715, 716, 718, 719, 720, 723, 725, 726, 727, 729, 730, 732, 733, 734, 735, 736, 746, 763, 765, 770, 771, 772, 773, 775, 776, 778, 780, 783, 791, 805, 810, 817, 821, 822, 823, 825, 827, 828, 829, 830, 832, 833, 834, 835, 836, 838, 842, 843, 844, 845, 846, 847, 848, 849, 850, 852, 853, 855, 856, 874, 877, 879, 880, 884, 886, 887, 889, 890, 895, 896, 897, 898, 904, 905, 906, 907, 908, 909, 910, 912, 913, 914, 915, 916, 917, 919, 920, 921, 922, 923, 924, 925, 926, 927, 928, 929, 930, 932, 935, 936, 938, 940, 941, 943, 945, 946, 949, 953, 956, 957, 963, 964, 970, 975, 980, 981, 985, 988, 991, 997, 998, 1002, 1004, 1008, 1009, 1011, 1012, 1013, 1018, 1020, 1021, 1022, 1029, 1033, 1034, 1044, 1045, 1048, 1057, 1074, 1094, 1095, 1102, 1135, 1144, 1154, 1160, 1167, 1180, 1225, 1228, 1231, 1240, None, None]
sector_assignments = [0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 2, 1, 0, 3, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 3, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 1, 3, 0, 3, 1, 0, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 1, 3, 1, 0, 0, 0, 3, 0, 3, 0, 0, 0, 0, 3, 0, 1, 3, 1, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 3, 1, 0, 0, 0, 0, 0, 3, 0, 3, 0, 0, 1, 0, 2, 0, 1, 0, 1, 0, 0, 3, 1, 0, 0, 1, 1, 0, 3, 1, 1, 3, 3, 0, 1, 2, 0, 0, 0, 3, 0, 0, 0, 0, 3, 3, 1, 1, 0, 3, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 3, 3, 1, 0, 1, 0, 3, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 1, 0, 3, 3, 0, 3, 0, 3, 0, 0, 0, 0, 0, 3, 0, 3, 1, 3, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 3, 0, 1, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 0, 3, 0, 1, 0, 3, 1, 3, 0, 0, 0, 0, 0, 3, 0, 0, 3, 0, 1, 0, 3, 0, 3, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 3, 0, 0, 1, 0, 3, 0, 3, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3]
# Define colors for each sector (0-3)
color_map = {
    0: "yellow",
    1: "red",
    2: "blue",
    3: "purple"
}

# Generate PyMOL commands
pml_commands = [
    "# PyMOL script generated for sector visualization",
    "# Delete any existing sector selections first",
    "delete sector_*",
    "",
    "# Create and color sector selections"
]

# Group residues by sector
from collections import defaultdict
sector_residues = defaultdict(list)
for res, sector in zip(sector_indices, sector_assignments):
    sector_residues[sector].append(res)

# Generate selection and color commands for each sector
for sector, residues in sector_residues.items():
    # Create selection (group residues in chunks to avoid command line limits)
    chunk_size = 100
    for i in range(0, len(residues), chunk_size):
        chunk = residues[i:i+chunk_size]
        resi_str = "+".join(map(str, chunk))
        pml_commands.append(f"select sector_{sector}_part{i//chunk_size}, resi {resi_str}")
    
    # Combine all parts
    parts = [f"sector_{sector}_part{i//chunk_size}" for i in range(0, len(residues), chunk_size)]
    pml_commands.append(f"select sector_{sector}, {' or '.join(parts)}")
    
    # Delete temporary parts
    for part in parts:
        pml_commands.append(f"delete {part}")
    
    # Color the sector
    pml_commands.append(f"color {color_map[sector]}, sector_{sector}")
    pml_commands.append(f"show sticks, sector_{sector}")
    pml_commands.append("")

# Add some visualization improvements
pml_commands.extend([
    "# Visualization settings",
    "set stick_radius, 0.2",
    "set sphere_scale, 0.3",
    "bg_color white",
    "set ray_opaque_background, off",
    ""
])

# Write to a PyMOL script file
with open("sectors.pml", "w") as f:
    f.write("\n".join(pml_commands))

print("PyMOL script 'sectors.pml' generated successfully!")