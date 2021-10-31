def find_all_L_sets():
    all_L_sets = []
    keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16', 'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
    for key1 in keys:
        incident_1 = set(get_incident_keys(key1))
        non_incident_1 = set(get_non_incident_keys(key1))
        for key2 in incident_1:
            incident_2 = set(get_incident_keys(key2))
            non_incident_2 = set(get_non_incident_keys(key2))
            for key3 in non_incident_1.intersection(incident_2):
                incident_3 = set(get_incident_keys(key3))
                non_incident_3 = set(get_non_incident_keys(key3))
                for key4 in incident_1.intersection(non_incident_2).intersection(incident_3):   
                    non_incident_4 = set(get_non_incident_keys(key4))
                    for key5 in non_incident_1.intersection(incident_2).intersection(non_incident_3).intersection(non_incident_4):                       
                        all_L_sets.append([key1, key2, key3, key4, key5])
    return all_L_sets