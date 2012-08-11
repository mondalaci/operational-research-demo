#!/usr/bin/env python
# -*- coding: iso-8859-2 -*-

import copy

import pygtk
pygtk.require('2.0')
import gobject
import gtk

# Application logic classes

class Matrix:
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.data = [[]] * x
        for x_i in range(x):
            self.data[x_i] = [None]*y

    def __repr__(self):
        s = ''

        for x_i in range(self.x):
            for y_i in range(self.y):
                s += `self.Get(x_i, y_i)` + ' '
            s += '\n'

        return s

    def Set(self, x, y, element):
        self.data[x][y] = element

    def Get(self, x, y):
        return self.data[x][y]

class SimplexMatrix(Matrix):
    ITER_OPTIMUM = 0
    ITER_NO_LIM = 1
    ITER_TRANSFORM = 2

    def __init__(self, matrix):
        self.x = x = matrix.x
        self.y = y = matrix.y

        Matrix.__init__(self, x, y)

        self.x_idx = range(1, x)
        self.y_idx = range(x, x+y-1)

        for i in range(x):
            for j in range(y):
                value = matrix.Get(i, j)
                self.Set(i, j, value)

    def __repr__(self):

        def p(s):
            return ('%5.3f' % s).rjust(8)

        def x(s):
            return ('x' + '%i' % s).rjust(8)

        s = ' ' * 8 + ' | '

        for i in range(self.y - 1):
            s += x(self.y_idx[i]) + ' '

        s += ' | \n'


        s += '-' * 9 + '+'
        s += '-' * ((9 * (self.y - 1)) + 2)
        s += '+' + '-' * 9 + '\n'


        for i in range(self.x - 1):
            s += x(self.x_idx[i]) + ' | '

            for j in range(self.y - 1):
                s += p(self.Get(i, j)) + ' '

            s += ' | ' + p(self.Get(i, self.y-1)) + '\n'


        s += '-' * 9 + '+'
        s += '-' * ((9 * (self.y - 1)) + 2)
        s += '+' + '-' * 9 + '\n'


        s += ' '*8 + ' | '

        for j in range(self.y - 1):
            s += p(self.Get(self.x - 1, j)) + ' '

        s += ' | ' + p(self.Get(self.x - 1, self.y-1)) + '\n'


        return s

    def Iterate(self):

        def get_cj_min():
            cj_pos = 0
            cj_old_min = self.Get(self.x-1, 0)

            for i in range(1, self.y-1):
                cj_new_val = self.Get(self.x-1, i)

                if cj_new_val < cj_old_min:
                    cj_old_min = cj_new_val
                    cj_pos = i

            if cj_old_min >= 0:
                return -1

            return cj_pos

        def get_lex_vect(row_idx, cj_idx):
            akj = self.Get(row_idx, cj_idx)
            vect = [self.Get(row_idx, self.y-1)/akj]

            row_dict = {}
            for row2_idx in self.x_idx:
                row_dict[row2_idx-1] = True

            for i in range(self.y-1 + self.x-1):
                if row_dict.has_key(i):

                    if self.x_idx[row_idx] == i+1:
                        val = 1/akj
                    else:
                        val = 0/akj

                else:
                    val = self.Get(row_idx-1, cj_idx)

                vect.append(val/akj)

            return vect

        def get_min_bk_akj_vectors(cj_idx):

            vectors = []
            akj = self.Get(0, cj_idx)
            bk_akj_min = 2**30

            if akj > 0:
                bk_akj_min = self.Get(0, self.x-1) / akj
                vect = (0, get_lex_vect(0, cj_idx))
                vectors.append(vect)

            for i in range(1, self.x-1):
                akj = self.Get(i, cj_idx)
                if akj > 0:
                    bk_akj = self.Get(i, self.x-1) / akj

                    if bk_akj < bk_akj_min:
                        bk_akj_min = bk_akj
                        vect = (i, get_lex_vect(i, cj_idx))
                        vectors = [vect]
                    elif bk_akj == bk_akj_min:
                        vect = (i, get_lex_vect(i, cj_idx))
                        vectors.append(vect)

            return vectors

        def get_min_vect(vectors):
            while len(vectors) > 1:
                for cur_idx in range(len(vectors[0][1])):
                    values = [vector[1][cur_idx] for vector in vectors]
                    min_value = min(values)
                    vectors = [vector for vector in vectors if vector[1][cur_idx]==min_value]

            return vectors[0]

        def do_simplex(col, row, self):
            smatrix = copy.deepcopy(self)

            akj = self.Get(row, col)

            val = self.x_idx[row]
            smatrix.y_idx[col] = val

            val = self.y_idx[col]
            smatrix.x_idx[row] = val

            val = 1 / akj
            smatrix.Set(row, col, val)

            for i in range(self.y):
                if i != row:
                    val = self.Get(row, i) / akj
                    smatrix.Set(row, i, val)

            for i in range(self.x):
                if i != col:
                    val = - self.Get(i, col) / akj
                    smatrix.Set(i, col, val)

            for i in range(self.x):
                for j in range(self.y):
                    if i != row and j != col:
                        val = self.Get(i, j) - self.Get(i, col)*self.Get(row, j)/akj
                        smatrix.Set(i, j, val)

            return smatrix


        cj_idx = get_cj_min()

        if cj_idx == -1:
            return (self.ITER_OPTIMUM, None, None)

        vectors = get_min_bk_akj_vectors(cj_idx)

        if len(vectors) == 0:
            return (self.ITER_NO_LIM, None, None)

        row = get_min_vect(vectors)[0]

        blurb = u'Iteráció, a generálóelem: %i, %i.\n\n' % (row+1, cj_idx+1)
        smatrix = do_simplex(cj_idx, row, self)

        return (self.ITER_TRANSFORM, smatrix, blurb)

    def SolveExercise(self):
        smatrix = self
        iter_state = self.ITER_TRANSFORM
        iter_num = 1
        text = ''

        while iter_state == self.ITER_TRANSFORM:
            text += u'%i. lépés:\n\n' % iter_num
            text += repr(smatrix) + '\n'
            iter_state, smatrix, blurb= smatrix.Iterate()

            if iter_state == self.ITER_TRANSFORM:
                text += blurb

            iter_num += 1

        if iter_state == self.ITER_OPTIMUM:
            text += u'A feladat megoldása optimális.'
        elif iter_state == self.ITER_NO_LIM:
            text += u'A feladat nemkorlátos a lehetséges megoldások halmazán.'

        return text

class Form:
    pass

class Exercise:
    def GetNameList(self, num, name):
        lst = []
        for i in range(1, num+1):
            lst.append(`i` + '. ' + name)
        return lst

    def GetSimplexMatrix(self, form):
        matrix = form.matrix

        for i in range(matrix.y-1):
            value = matrix.Get(matrix.x - 1, i)
            matrix.Set(matrix.x - 1, i, -value)

        return SimplexMatrix(matrix)


class Exercise1(Exercise):
    def GetMatrixWidget(self, form):
        form.ladies_num = form.input_fields[0]

        col_names = self.GetNameList(form.ladies_num, u'hölgy') + [u'Preferenciák']
        row_names = [u'Szépség (1-10)', u'Jellem (1-10)',
                     u'Megközelíthetőség (1-10)', u'Rászánt idő']
        matrix_widget = MatrixWidget(form.ladies_num+1, 4, row_names, col_names, True)

        return matrix_widget


class Exercise2(Exercise):
    def GetMatrixWidget(self, form):
        form.items_num = form.input_fields[0]
        form.metals_num = form.input_fields[1]

        row_names = self.GetNameList(form.items_num, u'Portéka') + \
                    [u'Rendelkezésre álló mennyiség']
        col_names = self.GetNameList(form.metals_num, u'Érc fajta') + \
                    [u'Egységnyi nyereség']
        matrix_widget = MatrixWidget(form.metals_num+1, form.items_num+1,
                                      row_names, col_names, True)

        return matrix_widget

class Exercise3(Exercise):
    def GetMatrixWidget(self, form):
        form.items_num = form.input_fields[0]
        form.metals_num = form.input_fields[1]

        row_names = self.GetNameList(form.items_num, u'Táplálék') + \
                    [u'Szükséges mennyiség']
        col_names = self.GetNameList(form.metals_num, u'Alapanyag') + \
                    [u'Egységnyi mennyiség ára']
        matrix_widget = MatrixWidget(form.metals_num+1, form.items_num+1,
                                      row_names, col_names, True)

        return matrix_widget

# Specialized widget classes

class Entries(gtk.Table):
    def __init__(self, content):
        length = len(content)
        gtk.Table.__init__(self, 2, length+1)
        self.set_row_spacings(4)
        self.input_widgets = []

        for row_i in range(length):
            row_name = content[row_i]

            prefix_label = gtk.Label(row_name + u':')
            prefix_label.set_padding(4, 0)
            prefix_label.set_alignment(1, -1)
            adjustment = gtk.Adjustment(1, 1, 65535, 1)
            spin_button = gtk.SpinButton(adjustment)
            self.input_widgets.append(spin_button)

            self.attach(prefix_label, 0, 1, row_i, row_i+1)
            self.attach(spin_button, 1, 2, row_i, row_i+1)

        button = gtk.Button(u'Tovább')
        button.connect('clicked', self.OnClicked)
        self.attach(button, 0, 3, row_i+1, row_i+2)

    def OnClicked(self, widget):
        form_data = Form()
        form_data.input_fields = [int(widget.get_value()) for widget in
                                  self.input_widgets]
        self.emit('submit_form', form_data)

class Border(gtk.VBox):
    def __init__(self, child):
        gtk.VBox.__init__(self)
        hbox = gtk.HBox()
        hbox.pack_start(child, False, False)
        hbox.pack_start(gtk.Label())

        self.pack_start(hbox, False, False)
        self.pack_start(gtk.Label())

class MatrixWidget(gtk.Table):
    def __init__(self, x, y, row_names=None, col_names=None, bottom_right_off=False):
        self.bottom_right_off = bottom_right_off
        gtk.Table.__init__(self, x+1, y+1)
        self.widget_matrix = Matrix(x, y)

        if row_names:
            for i in range(y):
                label = gtk.Label(row_names[i])
                label.set_padding(4,0)
                label.set_line_wrap(True)
                self.attach(label, i+1, i+2, 0, 1)

        if col_names:
            for i in range(x):
                label = gtk.Label(col_names[i])
                label.set_padding(4, 0)
                label.set_alignment(1, -1)
                label.set_line_wrap(True)
                self.attach(label, 0, 1, i+1, i+2, yoptions=gtk.EXPAND|gtk.FILL)

        for x_i in range(x):
            for y_i in range(y):
                entry = gtk.Entry()

                if x_i == x-1 and y_i == y-1 and bottom_right_off:
                    entry.set_sensitive(False)
                    entry.set_text(u'(nem érvényes)')

                entry.set_width_chars(8)
                self.widget_matrix.Set(x_i, y_i, entry)
                self.attach(entry, y_i+1, y_i+2, x_i+1, x_i+2)

    def GetMatrix(self):
        x = self.widget_matrix.x
        y = self.widget_matrix.y
        matrix = Matrix(x, y)

        try:
            for x_i in range(x):
                for y_i in range(y):
                    entry = self.widget_matrix.Get(x_i, y_i)

                    if x_i == x-1 and y_i == y-1 and self.bottom_right_off:
                        value = 0.0
                    else:
                        value = float(entry.get_text())

                    matrix.Set(x_i, y_i, value)
        except ValueError:
            return u'A (%i,%i) pozicióban szereplo érték nem szám!' % (x_i+1, y_i+1)

        return matrix

# Application Logic classes

class Application:
    def Init(self):
        gobject.signal_new('submit_form', Entries, gobject.SIGNAL_RUN_LAST,
                           gobject.TYPE_NONE, [gobject.TYPE_PYOBJECT])

        self.widget_stack = []

        self.InitGUI()

    def InitGUI(self):
        labels = [
            u'1. feladat:\nKaszanova',
            u'2. feladat:\nKovács',
            u'3. feladat:\nDiéta'
        ]

        hbox = gtk.HBox()

        for label_idx in range(len(labels)):
            button = gtk.Button(labels[label_idx])
            button.connect('clicked', self.OnExerciseStart, label_idx)
            hbox.pack_start(button)

        self.vbox = gtk.VBox()

        border = Border(self.vbox)

        sw = gtk.ScrolledWindow()
        sw.set_policy(gtk.POLICY_AUTOMATIC, gtk.POLICY_AUTOMATIC)
        sw.add_with_viewport(border)

        vbox = gtk.VBox()
        vbox.pack_start(hbox, False)
        vbox.pack_start(sw)

        window = gtk.Window()
        window.connect('destroy', self.OnQuit)
        window.set_title(u'Operációkutatás feladat')
        window.set_default_size(600, 500)
        window.add(vbox)
        window.show_all()

        gtk.main()

    def AddWidget(self, new_level, child):
        cur_level = len(self.widget_stack)
        add_level = cur_level - new_level

        for i in range(add_level):
            widget = self.widget_stack.pop()
            self.vbox.remove(widget)

        self.widget_stack.append(child)
        self.vbox.pack_start(child)
        self.vbox.show_all()

    def OnExerciseStart(self, widget, exercise_num):
        entries = [
            [u'Hölgyek száma'],
            [u'Portékák száma', u'Érc fajták száma'],
            [u'Táplálékok száma', u'Alapanyagok száma']
        ]

        self.exercise_num = exercise_num
        entries = Entries(entries[exercise_num])
        entries.connect('submit_form', self.OnSubmitForm1)

        border = Border(entries)
        self.AddWidget(0, border)

    def OnSubmitForm1(self, widget, form):
        self.AddMatrixWidget(form)

    def AddMatrixWidget(self, form):
        matrix_widget = exercises[self.exercise_num].GetMatrixWidget(form)

        vbox = gtk.VBox()
        vbox.pack_start(matrix_widget)
        button = gtk.Button(u'Eredmény')
        button.connect('clicked', self.OnSubmitForm2)
        vbox.pack_start(button)
        border = Border(vbox)
        self.AddWidget(1, border)
        self.form = form
        self.form.matrix_widget = matrix_widget

    def OnSubmitForm2(self, widget):
        matrix = self.form.matrix_widget.GetMatrix()
        self.form.matrix = matrix

        if type(matrix) == unicode:
            text = matrix
        else:
            simplex_matrix = exercises[self.exercise_num].GetSimplexMatrix(self.form)
            text = simplex_matrix.SolveExercise()

        self.AddTextView(text)

    def AddTextView(self, text):
        textbuffer = gtk.TextBuffer()
        textbuffer.insert_at_cursor(text)
        textbuffer.create_tag('fixed-width', font='MonoSpace 10')
        start_iter = textbuffer.get_start_iter()
        end_iter = textbuffer.get_end_iter()
        textbuffer.apply_tag_by_name('fixed-width', start_iter, end_iter)

        textview = gtk.TextView(textbuffer)
        textview.set_editable(False)
        textview.set_cursor_visible(False)

        self.AddWidget(2, textview)


    def OnQuit(self, widget):
        gtk.main_quit()


exercise1 = Exercise1()
exercise2 = Exercise2()
exercise3 = Exercise3()
exercises = [exercise1, exercise2, exercise3]
application = Application()

application.Init()
