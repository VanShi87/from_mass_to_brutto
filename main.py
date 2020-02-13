from tkinter import *

#source: http://moltensalt.org/references/static/downloads/pdf/stable-isotopes.pdf
exact_mass = {'C': 12.000000,
              'H': 1.007825,
              'N': 14.003074,
              'O': 15.994915,
              'S': 31.972071,
              'P': 30.973762,
              'F': 18.998403,
              'Cl': 34.968853,
              'Br': 80.916291,
              'I': 126.904468,
              'H+': 1.007825 - 0.000549,
              'Na+': 22.989770 - 0.000549,
              'K+': 39.963707 - 0.000549
              }
list_of_cations = ['H+', 'Na+', 'K+']
atom_list = []
Nmin = {}
Nmax = {}


def yield_number_of_atoms(atom, brutto, molecular_mass_rest):
    new_brutto = []
    new_molecular_mass_rest = []
    Nmax_theor = int(molecular_mass_rest / exact_mass[atom])
    if (Nmax[atom] is None) or (Nmax[atom] > Nmax_theor):
        Nmaximum = Nmax_theor
    else:
        Nmaximum = Nmax[atom]
    for i in range(Nmin[atom], Nmaximum + 1):
        if i == 0:
            new_brutto.append(brutto)
        else:
            new_brutto.append(brutto + atom + str(i))
        new_molecular_mass_rest.append(molecular_mass_rest - exact_mass[atom] * i)
    return new_brutto, new_molecular_mass_rest

def last_step(last_atom, brutto, molecular_mass_rest):
    for i, mass_rest in enumerate(molecular_mass_rest):
        N_last_atom = int(mass_rest / exact_mass[last_atom])
        if (Nmax[last_atom] is None and Nmin[last_atom]<=N_last_atom) or (Nmax[last_atom]>= N_last_atom and Nmin[last_atom]<=N_last_atom):
            brutto[i] += last_atom + str(N_last_atom)
            molecular_mass_rest[i] -= exact_mass[last_atom]*N_last_atom

def define_brutto():
    # инстанциируем список брутто-формул списком с единственным элементом - пустой строкой
    brutto = ['', ]
    molecular_mass = float(mol_mass_entry.get().replace(',', '.')) # replace - на случай использования запятой в качестве десятичного разделителя
    if unit_var.get():                                                                        # если погрешность выражена в миллионных долях,
        threshold = molecular_mass / 1000000 * float(threshold_entry.get().replace(',', '.')) # ее переводим в атомные единицы массы
    else:
        threshold = float(threshold_entry.get().replace(',', '.'))
    cation = list_of_cations[cation_var.get()]
    molecular_mass_rest = [molecular_mass + threshold - exact_mass[cation], ]
    element_list = [('C', C_var, C_min, C_max),
                 ('N', N_var, N_min, N_max),
                 ('O', O_var, O_min, O_max),
                 ('S', S_var, S_min, S_max),
                 ('P', P_var, P_min, P_max),
                 ('F', F_var, F_min, F_max),
                 ('Cl', Cl_var, Cl_min, Cl_max),
                 ('Br', Br_var, Br_min, Br_max),
                 ('I', I_var, I_min, I_max),
                 ('H', H_var, H_min, H_max),
                ]

    for element, element_var, element_min, element_max in element_list:
        if element_var.get():
            atom_list.append(element)
            try:
                Nmin[element] = int(element_min.get())
            except:
                Nmin[element] = 0
            try:
                Nmax[element] = int(element_max.get())
            except:
                Nmax[element] = None

    last_atom = atom_list.pop()

    for atom in atom_list:
        j = len(brutto)-1
        for i in range(j+1):
            brutto_new, molecular_mass_rest_new = yield_number_of_atoms(atom, brutto[j-i], molecular_mass_rest[j-i])
            del brutto[j-i]                   # эти значения либо дублируются (если Nmin[atom] == 0),
            del molecular_mass_rest[j-i]      # либо невозможны (если Nmin[atom] > 0)
            brutto.extend(brutto_new)
            molecular_mass_rest.extend(molecular_mass_rest_new)

    last_step(last_atom, brutto, molecular_mass_rest)



    for i, mass_rest in enumerate(molecular_mass_rest):
        if mass_rest <= 2 * threshold:
            print(brutto[i], mass_rest - threshold)

    atom_list.clear()  # обязательно (даже при реализации в форме множества)
    Nmin.clear()       # необязательно (неперезаписываемые значения в расчете не участвуют,
    Nmax.clear()       # однако могут помешать при тестировании)




window = Tk()
window.title('Расчет брутто-формулы на основе точного значения массы молекулы')
window.geometry('640x480')

mol_mass_label = Label(window, text='Молекулярная\nмасса, а.е.м.:')
mol_mass_label.grid(column=0, row=0)

threshold_label = Label(window, text='Допустимая\nпогрешность:')
threshold_label.grid(column=2, row=0)

mol_mass_entry = Entry(window, width=10)
mol_mass_entry.grid(column=1, row=0)

threshold_entry = Entry(window, width=10)
threshold_entry.grid(column=3, row=0)

unit_var = BooleanVar()
unit_var.set(0)
unit1 = Radiobutton(text='а.е.м.', variable=unit_var, value=0)
unit1.grid(column=2, row=1)
unit2 = Radiobutton(text='ppm', variable=unit_var, value=1)
unit2.grid(column=3, row=1)

atom_label = Label(window, text='Тип атома')
atom_label.grid(column=0, row=3)
min_label = Label(window, text='Мин. кол-во')
min_label.grid(column=1, row=3)
max_label = Label(window, text='Макс. кол-во')
max_label.grid(column=2, row=3)

C_var = BooleanVar()
C_var.set(0)
C_label = Checkbutton(window, text='C', variable=C_var, onvalue=1, offvalue=0)
C_label.grid(column=0, row=4)
C_min = Entry(window, width=10)
C_min.grid(column=1, row=4)
C_max = Entry(window, width=10)
C_max.grid(column=2, row=4)

N_var = BooleanVar()
N_var.set(0)
N_label = Checkbutton(window, text='N', variable=N_var, onvalue=1, offvalue=0)
N_label.grid(column=0, row=5)
N_min = Entry(window, width=10)
N_min.grid(column=1, row=5)
N_max = Entry(window, width=10)
N_max.grid(column=2, row=5)

O_var = BooleanVar()
O_var.set(0)
O_label = Checkbutton(window, text='O', variable=O_var, onvalue=1, offvalue=0)
O_label.grid(column=0, row=6)
O_min = Entry(window, width=10)
O_min.grid(column=1, row=6)
O_max = Entry(window, width=10)
O_max.grid(column=2, row=6)

S_var = BooleanVar()
S_var.set(0)
S_label = Checkbutton(window, text='S', variable=S_var, onvalue=1, offvalue=0)
S_label.grid(column=0, row=7)
S_min = Entry(window, width=10)
S_min.grid(column=1, row=7)
S_max = Entry(window, width=10)
S_max.grid(column=2, row=7)

P_var = BooleanVar()
P_var.set(0)
P_label = Checkbutton(window, text='P', variable=P_var, onvalue=1, offvalue=0)
P_label.grid(column=0, row=8)
P_min = Entry(window, width=10)
P_min.grid(column=1, row=8)
P_max = Entry(window, width=10)
P_max.grid(column=2, row=8)

F_var = BooleanVar()
F_var.set(0)
F_label = Checkbutton(window, text='F', variable=F_var, onvalue=1, offvalue=0)
F_label.grid(column=0, row=9)
F_min = Entry(window, width=10)
F_min.grid(column=1, row=9)
F_max = Entry(window, width=10)
F_max.grid(column=2, row=9)

Cl_var = BooleanVar()
Cl_var.set(0)
Cl_label = Checkbutton(window, text='Cl', variable=Cl_var, onvalue=1, offvalue=0)
Cl_label.grid(column=0, row=10)
Cl_min = Entry(window, width=10)
Cl_min.grid(column=1, row=10)
Cl_max = Entry(window, width=10)
Cl_max.grid(column=2, row=10)

Br_var = BooleanVar()
Br_var.set(0)
Br_label = Checkbutton(window, text='Br', variable=Br_var, onvalue=1, offvalue=0)
Br_label.grid(column=0, row=11)
Br_min = Entry(window, width=10)
Br_min.grid(column=1, row=11)
Br_max = Entry(window, width=10)
Br_max.grid(column=2, row=11)

I_var = BooleanVar()
I_var.set(0)
I_label = Checkbutton(window, text='I', variable=I_var, onvalue=1, offvalue=0)
I_label.grid(column=0, row=12)
I_min = Entry(window, width=10)
I_min.grid(column=1, row=12)
I_max = Entry(window, width=10)
I_max.grid(column=2, row=12)

H_var = BooleanVar()
H_var.set(0)
H_label = Checkbutton(window, text='H', variable=H_var, onvalue=1, offvalue=0)
H_label.grid(column=0, row=13)
H_min = Entry(window, width=10)
H_min.grid(column=1, row=13)
H_max = Entry(window, width=10)
H_max.grid(column=2, row=13)

mol_mass_label = Label(window, text='Катион:')
mol_mass_label.grid(column=5, row=0)
cation_var = IntVar()
cation_var.set(0)
cation1 = Radiobutton(text='H+', variable=cation_var, value=0)
cation1.grid(column=4, row=1)
cation2 = Radiobutton(text='Na+', variable=cation_var, value=1)
cation2.grid(column=5, row=1)
cation3 = Radiobutton(text='K+', variable=cation_var, value=2)
cation3.grid(column=6, row=1)

go_button = Button(window, text='Рассчитать', font=32, command=define_brutto)
go_button.grid(column=4, row=3, columnspan=3, rowspan=3)
window.mainloop()
print('end of the programm')



