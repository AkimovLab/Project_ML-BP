import os

def check(prefix,suffix,mini,maxi):
# The function creates the continuity regions of existing files
# prefix - prefix of the file to be found
# suffix - suffix of the file to be found
# mini - minimal index of the file to be found
# maxi - maximal index of the file to be found


    lst = []
    stop = 0
    i = mini
    while i<=maxi:
        file_path = prefix + str(i) + suffix       
        if not os.path.exists(file_path):
            lst.append(i)
        i = i + 1


    print("Missing files are = ", lst)

# Edit parameters of the function called below
# argument #1 - prefix of the file to be found
# argument #2 - suffix of the file to be found
# argument #3 - minimal index of the file to be found
# argument #4 - maximal index of the file to be found
print ("\n")
check(os.getcwd()+"/res/S_ks_","_re",5000,7000)
#check("res_mb_sp/Hvib_ci_", "_re", 5000, 7000)
print ("\nFinished Checking\n")
