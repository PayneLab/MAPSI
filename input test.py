print("Welcome to API_dev!\n")
input("How are you?")
input_folder_path = input("Input-file folder path: ")
if not input_folder_path.contains("\\"):
    print("Please enter path with \\")
    input_folder_path = input("Input-file folder path: ")