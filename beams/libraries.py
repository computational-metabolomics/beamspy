from collections import OrderedDict
import pandas

def _sort_nested_ordered_dict(od):
    res = OrderedDict()
    for k, v in sorted(od.items()):
        if isinstance(v, dict):
            res[k] = _sort_nested_ordered_dict(v)
        else:
            res[k] = v
    return res


class Atoms:
    def __init__(self):

        self.lib = OrderedDict()
        self.lib["C"] = 12.000000
        #self.lib["(13C)"] = 13.003355
        self.lib["H"] = 1.007825
        self.lib["N"] = 14.003074
        self.lib["O"] = 15.994915
        self.lib["P"] = 30.973763
        self.lib["S"] = 31.972072
        #self.lib["(34S)"] = 33.967868
        #self.lib["[NaCl]"] = 22.989770 + 34.968853
        #self.lib["[Na(37Cl)]"] = 22.989770 + 36.965903

        self.limits = OrderedDict()
        self.limits["C"] = range(0, 35)
        #self.limits["(13C)"] = range(0,2)
        self.limits["H"] = range(0, 73)
        self.limits["N"] = range(0, 16)
        self.limits["O"] = range(0, 20)
        self.limits["P"] = range(0, 8)
        self.limits["S"] = range(0, 9)
        #self.limits["(34S)"] = range(0,2)
        #self.limits["[NaCl]"] = range(0,2)
        #self.limits["[Na(37Cl)]"] = range(0,2)

        self.valence = {'C':4, 'H':1, 'N':3, 'O':2, 'P':3, 'S':2}

    def add(self, name, mass):
        self.lib[name] = mass
        self.lib = OrderedDict(sorted(self.lib.items(), key=lambda x:x[1]))

    def delete(self, name):
        if name == "*":
            for name in self.lib.keys():
                del self.lib[name]
        else:
            del self.lib[name]

    def __str__(self):
        out = "Atoms/Molecules in library\n"
        for key in self.lib:
            out += "%s\t%s\n" % (key, self.lib[key])
        out += "\n"
        out += "Atom limits\n"
        for key in self.limits:
            out += "%s\t%s-%s\n" % (key, min(self.limits[key]), max(self.limits[key]))
        return out

    def sort_length(self, items = []):

        def bylength(word1, word2):
            return len(word1) - len(word2)

        self.sort_atoms_length = self.lib.keys()
        self.sort_atoms_length.sort(cmp=bylength)
        return self.sort_atoms_length


class Adducts:
    def __init__(self, ion_mode=None, e=0.0005486):

        self.e = e
        if ion_mode == "pos":
            self.lib = OrderedDict()
        elif ion_mode == "neg":
            self.lib = OrderedDict()
        elif ion_mode is None:
            self.lib = OrderedDict()

    def add(self, name, mass):
        self.lib[name] = mass
        self.lib = OrderedDict(sorted(self.lib.items(), key=lambda x:x[1]))

    def delete(self, name):
        if name == "*":
            self.lib = OrderedDict()
        else:
            if name in self.lib:
                self.lib.remove(name)
            else:
                print "Entry not in library"

    def __str__(self):
        out = "Adducts in library\n"
        out += "-----------------\n"
        out += "name\tmass\n"
        for key in self.lib:
            out += "%s\t%s\n" % (key, self.lib[key])
        return out

    def sort_length(self, items = []):

        def bylength(word1, word2):
            return len(word1) - len(word2)

        self.sort_adducts_length = self.lib.keys()
        self.sort_adducts_length.sort(cmp=bylength)
        return self.sort_adducts_length


class Isotopes:

    def __init__(self, ion_mode=None):

        self.ion_mode = ion_mode
        self.lib = [OrderedDict([("C", {"abundance": 100.0}), ("(13C)", {"abundance": 1.1}), ("mass_difference", 1.003355)]),
                    OrderedDict([("S", {"abundance": 100.0}), ("(34S)", {"abundance": 4.21}), ("mass_difference", 1.995796)])]

        if self.ion_mode == "pos":
            self.lib.append(OrderedDict([("K", {"abundance": 100.0}), ("(41K)", {"abundance": 6.73}), ("mass_difference", 1.998117)]))
            #self.lib.append(OrderedDict([("(6Li)", {"abundance": 7.42}), ("Li", {"abundance": 1.0}), ("mass_difference", 1.000882)]))

        elif self.ion_mode == "neg":
            self.lib.append(OrderedDict([("Cl", {"abundance": 100.0}), ("(37Cl)", {"abundance": 24.23}), ("mass_difference", 1.997050)]))
        self.lib = sorted(self.lib, key=lambda k: k['mass_difference'])

    def add(self, label_x, label_y, mx_abundance, my_abundance, mass_difference):
        self.lib.append(OrderedDict([(label_x, {"abundance": float(mx_abundance)}),
                                     (label_y, {"abundance": float(my_abundance)}),
                                     ("mass_difference", mass_difference)]))
        self.lib = sorted(self.lib, key=lambda k: k['mass_difference'])

    def delete(self, label_x="*", label_y="*"):
        if label_x == "*" and label_y == "*":
            self.lib = []
        else:
            for item in self.lib:
                if label_x in item or label_y in item:
                    self.lib.remove(item)
            else:
                print "Entry not in library"

    def __str__(self):
        out = "Isotopes in library:\n"
        out += "--------------------------------------------\n"
        out += "label_x\tlabel_y\tmass_difference\tabundance_x\tabundance_y\n"
        for item in self.lib:
            label_x = item.items()[0][0]
            label_y = item.items()[1][0]
            out += "{}\t{}\t{}\t{}\t{}\n".format(label_x, label_y,
                                                     item["mass_difference"],
                                                     item[label_x]["abundance"], item[label_y]["abundance"])
        return out


class MultipleChargedDifferences:
    def __init__(self, ion_mode=None, e=0.0005486):

        self.e = e
        if ion_mode == "pos":
            self.lib = OrderedDict()
        elif ion_mode == "neg":
            self.lib = OrderedDict()
        elif ion_mode is None:
            self.lib = OrderedDict()

    def add(self, name, mass, charge, item_type):
        self.lib[name] = OrderedDict([("mass", float(mass)), ("charge", int(charge)), ("type", item_type)])
        self.lib = OrderedDict(sorted(self.lib.iteritems(), key=lambda x: x[1]['mass']))

    def delete(self, name):
        if name == "*":
            self.lib = OrderedDict()
        else:
            if name in self.lib:
                self.lib.remove(name)
            else:
                print "Entry not in library"

    def __str__(self):
        out = "Multiple charge ions in library\n"
        out += "-------------------------------\n"
        out += "name\tmass\tcharge\ttype\n"
        for key in self.lib:
            out += "{}\t{}\t{}\t{}\n".format(key, self.lib[key]["mass"], self.lib[key]["charge"], self.lib[key]["type"])
        return out

    def sort_length(self, items = []):

        def bylength(word1, word2):
            return len(word1) - len(word2)

        self.sort_adducts_length = self.lib.keys()
        self.sort_adducts_length.sort(cmp=bylength)
        return self.sort_adducts_length


class MultipleChargedIons:
    def __init__(self, ion_mode=None, e=0.0005486):

        self.e = e
        if ion_mode == "pos":
            self.lib = OrderedDict()
        elif ion_mode == "neg":
            self.lib = OrderedDict()
        elif ion_mode is None:
            self.lib = OrderedDict()

    def add(self, name, mass, charge):
        self.lib[name] = OrderedDict([("mass", float(mass)), ("charge", int(charge))])
        self.lib = OrderedDict(sorted(self.lib.iteritems(), key=lambda x: x[1]['mass']))

    def delete(self, name):
        if name == "*":
            self.lib = OrderedDict()
        else:
            if name in self.lib:
                self.lib.remove(name)
            else:
                print "Entry not in library"

    def __str__(self):
        out = "Multiple charge ions in library\n"
        out += "-------------------------------\n"
        out += "name\tmass\tcharge\n"
        for key in self.lib:
            out += "{}\t{}\t{}\n".format(key, self.lib[key]["mass"], self.lib[key]["charge"])
        return out

    def sort_length(self, items = []):

        def bylength(word1, word2):
            return len(word1) - len(word2)

        self.sort_adducts_length = self.lib.keys()
        self.sort_adducts_length.sort(cmp=bylength)
        return self.sort_adducts_length


class MassDifferences:

    def __init__(self, ion_mode=None):

        self.ion_mode = ion_mode
        self.lib = []

    def add(self, label_x, label_y, mass_difference, charge_x=1, charge_y=1):
        self.lib.append(OrderedDict([(label_x, {"charge": float(charge_x)}),
                                     (label_y, {"charge": float(charge_y)}),
                                     ("mass_difference", mass_difference)]))
        self.lib = sorted(self.lib, key=lambda k: k['mass_difference'])

    def delete(self, label_x="*", label_y="*"):
        if label_x == "*" and label_y == "*":
            self.lib = []
        else:
            for item in self.lib:
                if label_x in item or label_y in item:
                    self.lib.remove(item)
            else:
                print "Entry not in library"

    def __str__(self):
        out = "Mass differences in library:\n"
        out += "--------------------------------------------\n"
        out += "label_x\tlabel_y\tmass_difference\tcharge_x\tcharge_y\n"
        for item in self.lib:
            label_x, label_y = item.items()[0][0], item.items()[1][0]
            out += "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(label_x, label_y,
                                                         item[label_x]["charge"], item[label_y]["charge"],
                                                         item["mass_difference"])
        return out


class Formula:

    def __init__(self, formula, atoms, adducts, isotopes = []):

        self.atoms = atoms
        self.adducts = adducts
        self.isotopes = isotopes
        self.formula = formula

        self.count_atoms = {}
        self.count_adducts = {}

        self.exact_mass = 0.0
        self.mass = 0.0

        self.order_atoms = []

        self.formula_print_it = ""

    def count(self, debug=0):

        import re

        self.temp_formula = self.formula

        for atom in self.atoms.lib:
            self.count_atoms[atom] = 0
        for adduct in self.adducts.lib:
            self.count_adducts[adduct] = 0

        l_adducts = self.adducts.sort_length()
        l_adducts.reverse()

        for adduct in l_adducts:
            if adduct in self.temp_formula:
                self.temp_formula =  self.temp_formula.replace(adduct,"")
                self.count_adducts[adduct] += 1

        l_atoms = self.atoms.sort_length()
        l_atoms.reverse()

        for atom in l_atoms:
            str_pattern = []
            if "(" in atom or "[" in atom:

                num = len(atom)
                atom_search = atom.replace("(","\D").replace(")","\D").replace("[","\D").replace("]","\D")

                str_pattern.append("%s\d*" % (atom_search))
                pattern_char = re.compile(r'(%s)' % ("|".join(map(str,str_pattern))))
                result = pattern_char.findall(self.temp_formula)
                for item in result:
                    if atom in item:
                        self.temp_formula = self.temp_formula.replace(item, "")
                        count = item.replace(atom, "")
                        if len(count) == 0:
                            self.count_atoms[atom] += 1
                        else:
                            self.count_atoms[atom] += int(count)

            else:
                num = len(atom)
                str_pattern.append("%s\d*" % (atom))
                pattern_char = re.compile(r'(%s)' % ("|".join(map(str,str_pattern))))
                result = pattern_char.findall(self.temp_formula)
                for item in result:
                    if atom in item:
                        self.temp_formula = self.temp_formula.replace(item, "")
                        count = item.replace(atom, "")
                        if len(count) == 0:
                            self.count_atoms[atom] += 1
                        else:
                            self.count_atoms[atom] += int(count)

        if len(self.temp_formula) > 0:
            if debug == 1:
                print "Wrong Format Before:", self.formula, "After:", self.temp_formula
            for atom in self.atoms.lib:
                self.count_atoms[atom] = 0
            for adduct in self.adducts.lib:
                self.count_adducts[adduct] = 0
            return False

        return self.count_atoms, self.count_adducts

    def parent(self, isotopes = []):

        if isotopes != []:
            self.isotopes = isotopes

        self.count()

        if self.isotopes != []:
            for ip in self.isotopes.lib:
                if str(ip.items()[1][0]) in self.count_atoms:
                    self.count_atoms[ip[0]] += self.count_atoms[ip[1]]
                    self.count_atoms[ip[1]] = 0

        set_to_zero = []
        total_atoms = 0
        for a in self.count_atoms:
            if "[" and "]" in a and self.count_atoms[a] != 0:
                set_to_zero.append(a)
            else:
                total_atoms += self.count_atoms[a]
            if total_atoms > 0 and len(set_to_zero) > 0:
                for a in set_to_zero:
                    self.count_atoms[a] = 0

        return self.print_it(self.count_atoms)

    def calc_mass(self):

        self.count()

        self.exact_mass = 0.0
        self.mass = 0.0

        for atom in self.count_atoms:
            if self.count_atoms[atom] != 0:
                self.exact_mass += (self.count_atoms[atom] * (self.atoms.lib[atom]))

        self.mass += self.exact_mass
        for adduct in self.count_adducts:
            if self.count_adducts[adduct] != 0:
                self.mass += (self.count_adducts[adduct] * (self.adducts.lib[adduct]))

        return self.exact_mass, self.mass

    def HNOPS(self):

        self.parent() # update count_atoms

        self.HNOPS_rule = {"HC":0, "NOPSC":0}

        if (self.count_atoms['C']) == 0 or self.count_atoms['H'] == 0:
            self.HNOPS_rule["HC"]= 0
        elif (self.count_atoms['C']) == 0 and self.count_atoms['H'] == 0:
            self.HNOPS_rule["HC"]= 0
        elif (self.count_atoms['C']) != 0 and self.count_atoms['H'] != 0:
            if float(self.count_atoms['H'])/float((self.count_atoms['C'])) > 0 and float(self.count_atoms['H']/(self.count_atoms['C'])) < 6:
                self.HNOPS_rule["HC"]= 1
            if float(self.count_atoms['H'])/float((self.count_atoms['C'])) >= 6:
                self.HNOPS_rule["HC"]= 0

        #NOPS
        NOPS_check = []
        for atom in ['N','O','P','S']:
            try:
                check = float(float(self.count_atoms[atom]))/float((self.count_atoms['C']))
                NOPS_check.append(check)
            except ZeroDivisionError:
                NOPS_check.append(float(0))

        if NOPS_check[0] >= float(0) and NOPS_check[0] <= float(4) and NOPS_check[1] >= float(0) and NOPS_check[1] <= float(3) and NOPS_check[2] >= float(0) and NOPS_check[2] <= float(2) and NOPS_check[3] >= float(0) and NOPS_check[3] <= float(3):
            self.HNOPS_rule["NOPSC"] = 1
        if NOPS_check[0] > float(4) or NOPS_check[1] > float(3) or NOPS_check[2] > float(2) or NOPS_check[3] > float(3):
            self.HNOPS_rule["NOPSC"] = 0
        return self.HNOPS_rule

    def lewis_senior(self):

        self.parent() # update count_atoms

        self.lewis_senior_rule = {"lewis":0, "senior":0}
        #atoms_sum = sum(self.count_atoms.values())

        atoms_sum = sum(self.count_atoms.values())

        lewis_sum = 0

        for atom in self.atoms.valence:
            if atom == 'S':
                lewis_sum += self.atoms.valence[atom] * (self.count_atoms[atom])

            elif atom == 'C':
                lewis_sum += self.atoms.valence[atom] * (self.count_atoms[atom])

            else:
                lewis_sum += self.atoms.valence[atom] * self.count_atoms[atom]

        if lewis_sum%2 == 0:
            self.lewis_senior_rule["lewis"] = 1
        if lewis_sum%2 != 0:
            self.lewis_senior_rule["lewis"] = 0
        if lewis_sum >= ((atoms_sum-1)*2):
            self.lewis_senior_rule["senior"] = 1
        if lewis_sum < ((atoms_sum-1)*2):
            self.lewis_senior_rule["senior"] = 0
        self.count()
        atoms_sum = sum(self.count_atoms.values())
        return self.lewis_senior_rule

    def print_it(self, atom_count={}, adduct_count={}):

        if len(atom_count) > 0:
            self.count_atoms = atom_count
        if len(adduct_count) > 0:
            self.count_adducts = adduct_count
        if atom_count == {} and adduct_count == {}:
            self.count()

        self.formula_print_it = ""
        self.order_atoms = ["C", "(13C)", "H", "N", "O", "P", "S", "(34S)"]
        for atom in self.atoms.sort_length():
            if atom not in self.order_atoms:
                self.order_atoms.append(atom)

        for atom in self.order_atoms:
            try:
                if abs(self.count_atoms[atom]) > 1:
                    self.formula_print_it += atom + str(self.count_atoms[atom])
                if abs(self.count_atoms[atom]) == 1:
                    self.formula_print_it += atom
            except:
                pass

        for adduct in self.adducts.sort_length():
            try:
                if self.count_adducts[adduct] > 1:
                    self.formula_print_it += adduct + str(self.count_adducts[adduct])
                if self.count_adducts[adduct] == 1:
                    self.formula_print_it += adduct
            except:
                pass
        if self.formula_print_it == "":
            return False
        else:
            return self.formula_print_it

    def diff(self, inp_formula):

        self.count()

        self.inp_formula = Formula(inp_formula, self.atoms, self.adducts, self.isotopes)
        self.inp_formula.count()

        self.formulae_diff = ""

        for item in self.count_atoms.keys():
             self.inp_formula.count_atoms[item] -= self.count_atoms[item]

        for atom in self.atoms.sort_length():
            if self.inp_formula.count_atoms[atom] == 1:
                self.formulae_diff += "-" + str(atom)

            if self.inp_formula.count_atoms[atom] == -1:
                self.formulae_diff += str(atom)

            if self.inp_formula.count_atoms[atom] < -1:
                self.formulae_diff += str(atom) + str(abs(self.inp_formula.count_atoms[atom]))

            if self.inp_formula.count_atoms[atom] > 1:
                self.formulae_diff += "-" + str(atom) + str(self.inp_formula.count_atoms[atom])

        return self.formulae_diff
