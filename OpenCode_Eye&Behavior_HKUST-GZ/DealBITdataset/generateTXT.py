from backupFunctions import *

# Get a sorted list of directories in the current folder with numeric names
# Open, read, and write the content of the current text file to the output file
dir_list = sorted([f for f in os.listdir('.') if os.path.isdir(f) and f.isdigit()])

output_file = 'etgdata.txt'

for dir_name in dir_list:
    etgtxt_path = os.path.join(dir_name,'gazepoint')
    output_etgtxt = os.path.join(dir_name,output_file)
    etgtxtframes = sorted([f for f in os.listdir(etgtxt_path) if f.endswith('.txt')])

    with open(output_etgtxt, 'w') as etg_data:
        for txt in etgtxtframes:
            txt_path = os.path.join(etgtxt_path, txt)
            print(txt_path)
            with open(txt_path, 'r') as txt_data:
                data = txt_data.read()
                etg_data.write(data)