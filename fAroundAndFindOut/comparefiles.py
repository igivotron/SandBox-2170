
file1 = "output.obj"
file2 = "outputbis.obj"

with open(file1, 'r') as f1, open(file2, 'r') as f2:
    lines1 = f1.readlines()
    lines2 = f2.readlines()

differences = []
for i, (line1, line2) in enumerate(zip(lines1, lines2)):
    if line1 != line2:
        differences.append((i + 1, line1.strip(), line2.strip()))
if differences:
    print("Differences found:")
    for line_num, line1, line2 in differences:
        print(f"Line {line_num}:\n  File1: {line1}\n  File2: {line2}")

else:
    print("The files are identical.")
    