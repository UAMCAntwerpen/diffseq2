"""
(c) 2016 Joao Rodrigues
Adapted by Olivier Beyens (diffseq --> diffseq2)
"""

from __future__ import print_function
from re import A
from pymol import cmd


# Autocomplete
def _autocomplete():
    return cmd.Shortcut(cmd.get_names())
cmd.auto_arg[0]['diffseq'] = [_autocomplete, 'object name', '']
cmd.auto_arg[1]['diffseq'] = [_autocomplete, 'object name', '']


# Helper focunction : Object to all resids (contains duplicates)
def _object_to_all_resids(object):
    """
    DESCRIPTION

    Takes an object and returns a list of all the resids
    
    
    
    """
    
    _s = {'_r': []}
    cmd.iterate(f'{object}',
                    '_r.append(resi)', space=_s)
    all_resids_with_duplicates = _s['_r']

    all_resids = []
    seen = set()
    for item in all_resids_with_duplicates:
        if item in seen:
            pass
        else:
            all_resids.append(item)
            seen.add(item)
    return all_resids

def identify_blocks(l):
    """
    DESCRIPTION
    Divides up a list in blocks. Each entry of sequential numbers is transformed to a block. Example: [1,2,3,11,12,13] --> [[1,2,3],[11,12,13]]

    INPUT:
    l: list of numbers
    
    OUTPUT:
    a list of lists. Each inner list is a consecutive block.
    
    """
    index = 0
    output = []
    current_block = []
    while index != len(l):

        current_entry = l[index]
        int_current_entry = int(current_entry)
        if index == 0:
            current_block.append(l[index])
        else:
            
            
            if str(int_current_entry-1) == l[index-1]:
                current_block.append(current_entry)
            else:
                output.append(current_block)
                current_block = [current_entry]
                #special case for the last one
            if index == len(l)-1:
                output.append(current_block)
        index +=1


    return output

def search_color(object):
        _s = {'_r': []}
        cmd.iterate(f'{object} and name CA',
                '_r.append(color)', space=_s)
        return _s['_r'][0]


# Map id to resid
def couple_to_resid(couple,object1, object2):
    """
    DESCRIPTION
    Takes an aligned couple [(object,atomid),(object,atomid)] and convers it to the corresponding residue'id's. outputs residue id of object1 and object2, in that order.
    
    
    INPUT
    couple: list of two tuples. two tuples contain the object name and the atom id
    object1: string with object 1 name
    object2: string with object 2 name

    OUTPUT
    residue id of object 1 and object 2, in that order
    
    
    """

    entry1, entry2 = couple
    if entry1[0] == object1:
        objA, at_id_A = entry1
        objB, at_id_B = entry2
    elif entry1[0] == object2:
        objA, at_id_A = entry2
        objB, at_id_B = entry1
    else:
        print('ERRORCODE1: object names not found')       
    
    
    def _obj_and_id_to_resid(obj,at_id):
        _s = {'_r': []}
        cmd.iterate(f'{obj} and id {at_id}',
                '_r.append(resi)', space=_s)
        return _s['_r'][0]
    residA = _obj_and_id_to_resid(objA,at_id_A)
    residB = _obj_and_id_to_resid(objB,at_id_B)
    return residA, residB



# Map id to resn
def couple_to_resn(couple, object1,object2):
    """
    DESCRIPTION
    Takes an aligned couple [(object,atomid),(object,atomid)] and convers it to the corresponding residue names. outputs residue name of object1 and object2, in that order.
    
    
    INPUT
    couple: list of two tuples. two tuples contain the object name and the atom id
    object1: string with object 1 name
    object2: string with object 2 name

    OUTPUT
    residue name of object 1 and object 2, in that order
    
    
    """
    entry1, entry2 = couple
    if entry1[0] == object1:
        objA, at_id_A = entry1
        objB, at_id_B = entry2
    elif entry1[0] == object2:
        objA, at_id_A = entry2
        objB, at_id_B = entry1
    else:
        print('ERRORCODE1')
    
    def _obj_and_id_to_resn(obj,at_id):
        _s = {'_r': []}
        cmd.iterate(f'{obj} and id {at_id}',
                '_r.append(resn)', space=_s)
        return _s['_r'][0]
    resnA = _obj_and_id_to_resn(objA,at_id_A)
    resnB = _obj_and_id_to_resn(objB,at_id_B)
    return resnA,resnB

# Go from resid to residue name
def resid_to_resn(resid,obj):
    """
    DESCRIPTION
    Takes a resid of an object and returns the resname of the corresponding resid and resname.

    INPUT:
    resid: string of the resid
    obj: pymol object

    OUTPUT
    corresponding resname
    
    
    """

    _s = {'_r': []}
    cmd.iterate(f'{obj} and resid {resid}',
                '_r.append(resn)', space=_s)
    return _s['_r'][0]    

# Compare if two blocks are the same
def compare_two_blocks_if_same(block1,block2,object1,object2):
    """
DESCRIPTION
    This function compares two blocks of resids of EQUAL LENGTH to check if the resname sequence of these blocks are the same.

INPUT:
    block1: resids of object1
    block2: resids of object2
    object1: first pymol object
    object2: second pymol object

OUTPUT:
    True/False depending on wether the sequence is actually the same or not.
    
    """

    output = True
    index = 0
    while index != len(block1):
        resn_A = resid_to_resn(block1[index],object1)
        resn_B = resid_to_resn(block2[index],object2)

        if resn_A != resn_B:
            output = False
        else:
            pass
    
        index+=1

    return output


def diffseq2(objectA, objectB,howmuchdiff= 'maxdiff'):
    """
DESCRIPTION

    Highlight sequence differences between two proteins. 
    
    maxdiff/lessdiff controls what to do with parts of the sequence that are not considered aligned by the pymol alignment algorithm
    maxdiff: resiudes that are not considered aligned by pymol, but have the same residue name are considered as different. (default)
    lessdiff: resiudes that are not considered aligned by pymol, but have the same residue name can be considered as sufficiently different.
        In the case of lessdif we identify three different secenarios:
        1) the non aligned block has the same length in A and  in B and has the exact same residue names
        2) the non aligned block has the same length in A and  in B and has some different residue names 
        3) the non aligned block has  not the same length in A and  in B

        We have chosen to implement 'lessdiff' in such a way that only case 1 will not be included in the difference objects. Case 2 and 3 are seen as suffieciently different
USAGE

    diffseq2 protA, protB
    diffseq2 protA, protB, maxdiff 
    diffseq2 protA, pritB, lessdiff
    """

    # Step 0: Check the input if it is valid
    if howmuchdiff != 'maxdiff' and howmuchdiff != 'lessdiff':
        print('Please enter a correct value for maxdiff/lessdiff')

    #############################
    # Step 1: Align both proteins
    cmd.align(objectA, objectB, object=f'aln_{objectA}_{objectB}')
    raw_aln = cmd.get_raw_alignment(f'aln_{objectA}_{objectB}')

    # Initialise the list of resids that are in alignment
    A_resid_in_align = []
    B_resid_in_align = []
    matches_of_A = {} # keys: resids of A, entries resids of B

    #Initialise the list of differences to be visiualised
    A_resid_difference_list = []
    B_resid_difference_list = []

    #############################
    # Step 2: Compare the aligned residues to see if thezy are the same or not
    for aligned_couple in raw_aln:
        
        
        resn_A,resn_B = couple_to_resn(aligned_couple,objectA,objectB)
        resid_A,resid_B = couple_to_resid(aligned_couple,objectA,objectB)

        if resid_A not in A_resid_in_align:
            A_resid_in_align.append(resid_A)
        
        if resid_B not in B_resid_in_align:
            B_resid_in_align.append(resid_B)

        if resn_A != resn_B:


            if resid_A not in A_resid_difference_list:
                A_resid_difference_list.append(resid_A)
            
            if resid_B not in B_resid_difference_list:
                B_resid_difference_list.append(resid_B)
        
        # Also if residue names are not the same we must add them in the match list to use to identfy the residue before the non aligned block
        if resid_A not in matches_of_A:
            matches_of_A[resid_A] = resid_B


    #cmd.delete(f'aln_{objectA}_{objectB}')

    ############################# 
    # Step 3: Look for Diffetences in non-aligned residues




    # Step 3.1 look for all resids per protein
    resids_A = _object_to_all_resids(objectA)
    resids_B = _object_to_all_resids(objectB)



    # Step 3.2: Case A: all non aligned residues are to be considered different    
    if howmuchdiff == 'maxdiff':
        for res in resids_A:
            if res not in A_resid_in_align:
                A_resid_difference_list.append(res)

        for res in resids_B:
            if res not in B_resid_in_align:
                B_resid_difference_list.append(res)
    
    # Step 3.2: Case B: if the non aligned sequence is exactly the same, not considered different 
    elif howmuchdiff == 'lessdiff':
        
        # 3.2B.1 Compose a list of the non-aligned residues
        special_listA = []
        special_listB = []
        for res in resids_A:
            if res not in A_resid_in_align:
                special_listA.append(res)

        for res in resids_B:
            if res not in B_resid_in_align:
                special_listB.append(res)

        # 3.2B.2  we initialize the lists we already considered
        considered_blocklistB = []

        # 3.2B.3  we divide the non aligned parts in blocks:
        non_aligned_blocks_A =  identify_blocks(special_listA)
        non_aligned_blocks_B =  identify_blocks(special_listB)
        

        print(non_aligned_blocks_A)

        # 3.2B.4 for each block of A, we check if case 1 is valid
        for blockofA in non_aligned_blocks_A:
            #  see if there is a corresponding B block:
            # to do this: find the previous resid before the block starts, and check what the parner is in B
            # If hte partner in B is the residue before the start of the B block, we have found a matching block.
            first_resid_block_A = blockofA[0]
            first_resid_block_A_int = int(first_resid_block_A)
            previous_resid_of_A= str(first_resid_block_A_int-1)


            # Special case: we are dealing with the first non aligned block and this is at the start of both proteins

            ###### Addendum ########
            # There are fout possibilities for a the first non aligned block
            # 1. Both protein start with a non aligned block
            # 2. Protein A starts with a non alingned sequence, but not B
            # 3. Protein B starts with a non alingned sequence, but not A
            # 4. Both proteins start with an aligned block

            # In case 4, we should move to the general case
            # In case 3, we can also safely go to the general case. The first not aligned B block will not end up in the considered blocks, so it will end in the selection.
            # in cases 2 and 4: these cannot be handled by the general case, since looking for the previous match will fail and we cannot identify a previous resid.

            #Thus we need to find an if statement that includes cases 2 and 4: that is the case if protein A starts with a non alined block
            ######################
            if blockofA == non_aligned_blocks_A[0] and first_resid_block_A == resids_A[0] :
                ####### Addendum ######
                # Now we have to discriminate between case 2 and 4. In case of case 2, we should add the first non aligned block to the differences in any case.
                
                # In case 4, we will have to first check if both non aligned parts have the same length or not.
                # 
                ######################

                if non_aligned_blocks_B[0][0] == resids_B[0]:
                    if len(non_aligned_blocks_B[0]) == len(blockofA):
                        
                        comaprison = compare_two_blocks_if_same(blockofA,non_aligned_blocks_B[0], objectA,objectB)

                        if comaprison == True:
                            pass

                        else:
                            for x in blockofA:
                                A_resid_difference_list.append(x)
                        
                            for x in non_aligned_blocks_B[0]:
                                B_resid_difference_list.append(x)
                    

                    # not the same length: considered different enough:
                    else:
                        for x in blockofA:
                            A_resid_difference_list.append(x)
                        
                        for x in non_aligned_blocks_B[0]:
                            B_resid_difference_list.append(x)

                    # We have now properly considered the first block of B:
                    considered_blocklistB.append(non_aligned_blocks_B[0])

                else:
                    A_resid_difference_list.append(blockofA)

            # Otherwise: general case:
            else:
                previous_resid_of_B = matches_of_A[previous_resid_of_A]
                possible_first_resid_of_corresponding_B_block = str(int(previous_resid_of_B)+1)

                corresponding_b_block = False
                for blockofB in non_aligned_blocks_B:
                    if possible_first_resid_of_corresponding_B_block in blockofB:
                        corresponding_b_block = blockofB
            
                    else:
                        pass
            
                # if there is no corresponding B block, then we add theresidues of A block to diffrences:
                if corresponding_b_block == False:
                    for entry in blockofA:
                        A_resid_difference_list.append(entry)
            
                # When there is a corresponding block of B
                else:
                    if len(blockofA) == len(corresponding_b_block):
                        # Now check wether sequence is the same
                        comaprison = compare_two_blocks_if_same(blockofA,corresponding_b_block,objectA,objectB)

                        # Case 1
                        if comaprison == True:
                            pass

                        # Case 2
                        else:
                            for entry in blockofA:
                                A_resid_difference_list.append(entry)
                            for entry in corresponding_b_block:
                                B_resid_difference_list.append(entry)
                    
                    
                    
                    

                    # Then case 3 if not the same length:
                    else:
                        for entry in blockofA:
                            A_resid_difference_list.append(entry)
                        for entry in corresponding_b_block:
                            B_resid_difference_list.append(entry)
                
                    # Add that this B_block has been considered
                    considered_blocklistB.append(corresponding_b_block)

        for blockB in non_aligned_blocks_B:
            
            # check if the block has already been considered
            # if not, there cannot be a corresponding A block of the same length:
            if blockB not in considered_blocklistB:
                for entry in blockB:
                    B_resid_difference_list.append(entry)


            

    else:
        print('ERRORCODE 2: no valid howmuchdiff')


    #Convert resid difference lsit to a selection string for pymol
    def build_selection_string(difference_list):
        selection_string=''
        for entry in difference_list:
            selection_string += ',' + str(entry)
        selection_string = selection_string[1:] # remove first comma

        return selection_string
    
    selection_string_A = build_selection_string(A_resid_difference_list)
    selection_string_B = build_selection_string(B_resid_difference_list)

    cmd.select(f'{objectA}_diff_from_{objectB}',f'{objectA} and resid {selection_string_A}') # selection name, selection string
    
    cmd.select(f'{objectB}_diff_from_{objectA}',f'{objectB} and resid {selection_string_B}') # selection name, selection string


    cmd.create(f'{objectA}_vs_{objectB}',f'{objectA}_diff_from_{objectB}')
    cmd.create(f'{objectB}_vs_{objectA}',f'{objectB}_diff_from_{objectA}')

    cmd.do(f'hide everything, {objectA}_vs_{objectB}')
    cmd.do(f'hide everything, {objectB}_vs_{objectA}')

    cmd.show('surface',f'{objectA}_vs_{objectB}')
    cmd.show('surface',f'{objectB}_vs_{objectA}')

    cmd.set('transparency', '0.55')


    colorA = search_color(objectA)
    colorB = search_color(objectB)

    cmd.do(f'color {colorA},{objectA}_vs_{objectB}')
    cmd.do(f'color {colorB},{objectB}_vs_{objectA}')

    cmd.refresh()

cmd.extend('diffseq2', diffseq2)
