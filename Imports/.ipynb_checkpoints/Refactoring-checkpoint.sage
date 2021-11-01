#Tofix    
def find_all_permutations():
    all_perm = []
    keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16', 'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
    for key1 in keys:
        non_in_1 = set(get_non_incident_keys(key1))
        in_1 = set(get_incident_keys(key1))
        for key2 in non_in_1:
            non_in_2 = set(get_non_incident_keys(key2))
            in_2 = set(get_incident_keys(key2))
            for key3 in non_in_1.intersection(non_in_2):
                non_in_3 = set(get_non_incident_keys(key3))
                in_3 = set(get_incident_keys(key3))
                for key4 in non_in_1.intersection(non_in_2).intersection(non_in_3):
                    non_in_4 = set(get_non_incident_keys(key4))
                    in_4 = set(get_incident_keys(key4))
                    for key5 in non_in_1.intersection(non_in_2).intersection(non_in_3).intersection(non_in_4):
                        non_in_5 = set(get_non_incident_keys(key5))
                        in_5 = set(get_incident_keys(key5))
                        for key6 in non_in_1.intersection(non_in_2).intersection(non_in_3).intersection(non_in_4).intersection(non_in_5):
                            non_in_6 = set(get_non_incident_keys(key6))
                            in_6 = set(get_incident_keys(key6))
                            E = [key1, key2, key3, key4, key5, key6]
                            non_incident = [non_in_1, non_in_2, non_in_3, non_in_4, non_in_5, non_in_6]
                            incident = [in_1, in_2, in_3, in_4, in_5, in_6]
                            G = []
                            for i in range(6):
                                s = set([key for key in keys if key not in E])
                                for j in range(6):
                                    if j!=i:
                                        s.intersection_update(incident[j])
                                    else:
                                        s.intersection_update(non_incident[j])
                                G.append(list(s)[0])
                            F=[]
                            for i in range(5):
                                for j in range(i+1, 6):
                                    s = set([key for key in keys if key not in E and key not in G])
                                    for k in range(6):
                                        if i == k or j == k:
                                            s.intersection_update(incident[k])
                                        else:
                                            s.intersection_update(non_incident[k])
                                    F.append(list(s)[0])
                            all_perm.append(E+G+F)                           
    return all_perm