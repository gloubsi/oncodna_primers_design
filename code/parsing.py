from target import Target
import csv


def get_data_from_csv(name_of_file):
    """
    Parse the file containing all the information

    :return: a list of target objects
    """
    res = []
    with open(name_of_file, 'r') as f:
        reader = csv.reader(f)
        # skip the first line (header)
        next(reader)
        for row in reader:
            data = row[0].split()
            if data[11] == "Damaging":
                res.append(Target(data))
            elif data[11] == "Probably":
                data[11] = data[11] + " " + data[12]
                del data[12]
                res.append(Target(data))
            elif data[11] == "Potentially":
                data[11] = data[11] + " " + data[12]
                del data[12]
                res.append(Target(data))
            elif data[11] == "Unknown":
                res.append(Target(data))
    return res
