def change_line_with_number(line_number_start, line_number_end, number, file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for i in range(1, len(lines) + 1):
        if i >= line_number_start and i <= line_number_end:
            print(lines[i - 1])
            lines[i - 1] = str(number) + "\n"
            print(f"Line {i} in {file_path} has been changed to {number}.")
    
    with open(file_path, 'w') as file: 
        file.writelines(lines)

#usage:
line_number_start = int(input("Enter the line number to start: "))
line_number_end = int(input("Enter the line number to end: "))
number = 3000000000#int(input("Enter the number: "))
file_path = '/home/clemdoe/OpenFOAM/GeN-Foam/PhaseII_Ex2_04_6770_changepower/0/fluidRegion/powerDensity.fixedPower'

change_line_with_number(line_number_start, line_number_end, number, file_path)
