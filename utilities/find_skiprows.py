"""
Define functions to help read text files
"""

def find_skiprows(full_path):
    '''
    Find the first row where the text file starts with a number
    This allows the value of skiprows in text reading code to be 
    found automatically.


    WARNING: This function fails if the first number is a negative. The "-" sign doesn't convert to a float. 

    Parameters
    ----------
    full_path : STR
        Full path to file including file name and extension.

    Returns
    -------
    i : INT
        Number of rows until number is found.
    '''

    txt_file = open("%s" % (full_path),'rt')
    rows     = txt_file.readlines()
    
    test = False
    i    = 0
    while test == False:
        try:
            float(rows[i][0]) # Check if first character can be converted to float
            test = True
        except:
            i = i+1
    return i