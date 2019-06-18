import os
import shutil

def magicRename(folder_to_duplicate):
    # save the current path
    curPath = os.getcwd()
    # make a copy of the folder before processing the files
    all_dirs_in_path = [d for d in os.listdir(curPath) if os.path.isdir(d)]

    all_converted_dirs_in_path = [d for d in all_dirs_in_path if 'converted' in d]

    largest_index = max([0] + [int(d.split('_')[-1]) for d in all_converted_dirs_in_path])
    new_index = largest_index + 1
    new_name = folder_to_duplicate + '_converted_%d' %new_index
    shutil.copytree(os.path.join(curPath, folder_to_duplicate), os.path.join(curPath, new_name))

    filesList = os.listdir(os.path.join(curPath, new_name))
    os.chdir(os.path.join(curPath, new_name))
    for fileName in filesList:
        print "The file: [" + fileName + "] renamed to: [" + fileName.translate(None, "0123456789") + "]"
        os.rename(fileName, fileName.translate(None, "0123456789"))
        print "check your files now!"
    # back to old path
    os.chdir(curPath)

if __name__ == '__main__':
    magicRename(r'C:\Users\Zhongming\Dropbox\The D paper\data\D_simulations\0')