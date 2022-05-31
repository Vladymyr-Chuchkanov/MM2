import copy
from tkinter import ttk
from tkinter import *
from tkinter import messagebox
import numpy as np
import re
import sympy as sp
from scipy import integrate
import math
MathWords = 'sin cos arccos arcsin tg arctg arcch arcsh sh ch th log e Pi t pow'
def parse0(str0):
    str0 = str0.lower()
    calc0 = len(re.findall('\(', str0))
    if calc0 != len(re.findall('\)', str0)):
        return False;
    calc0 = re.findall('[a-zA-Z]+', str0)
    for el in calc0:
        print(re.search(el, MathWords))
        if re.search(el, MathWords) is None:
            return False
    str0 = re.sub(' ', '', str0)
    calc0 = re.sub('[a-zA-Z0-9+/^*,.\|\-\(\)]', '', str0)
    if len(calc0) > 0:
        return False
    return True
odn = 0
def parseEval(t):
    global flag
    global str000
    res = 0

    cop = re.sub('arctg', '&0', str000)
    cop = re.sub('math.', '', cop)
    cop = re.sub('tg', '&1', cop)
    cop = re.sub('th', '&2', cop)
    cop = re.sub('t', str(t), cop)
    cop = re.sub('&0', 'math.atan', cop)
    cop = re.sub('&1', 'math.tan', cop)
    cop = re.sub('&2', 'math.tanh', cop)
    cop = re.sub('sin', 'math.sin', cop)
    cop = re.sub('cos', 'math.cos', cop)
    cop = re.sub('e', 'math.e', cop)
    cop = re.sub('Pi', 'math.pi', cop)
    cop = re.sub('arcsin', 'math.asin', cop)
    cop = re.sub('arccos', 'math.acos', cop)
    cop = re.sub('arcch', 'math.acosh', cop)
    cop = re.sub('arcsh', 'math.asinh', cop)
    cop = re.sub('lg', 'math.log', cop)
    cop = re.sub('pow', 'math.pow', cop)
    try:
        res = eval(cop)
    except:
        messagebox.showinfo("Error", "wrong function "+cop+" !")
        return 0
    return res

def matrixMulti(A,B):
    resM = len(B[0])
    resN = len(A)
    a2 = len(B)
    a1 = len(A[0])
    if a1!=a2:
        messagebox.showinfo("Error", "wrong size of Matrix to multi!")
        return 0
    res =[]
    for i in range(resN):
        res.append([])
        for j in range(resM):
            temp =''
            for a in range(a1):
                if a!=0:
                    temp+='+'
                temp+='(('+str(A[i][a])+')*('+str(B[a][j])+'))'
            res[i].append(temp)
    return res
def matrixAdd(A,B):
    resM = len(B[0])
    resN = len(A)
    a2 = len(B)
    a1 = len(A[0])
    if a1!=resM and a2!=resN:
        messagebox.showinfo("Error", "wrong size of Matrix to add!")
        return 0
    res = []
    for i in range(resN):
        res.append([])
        for j in range(resM):
            temp = ''
            temp += '((' + str(A[i][j]) + ')+(' + str(B[i][j]) + '))'
            res[i].append(temp)
    return res
def matrixMinus(A,B):
    resM = len(B[0])
    resN = len(A)
    a2 = len(B)
    a1 = len(A[0])
    if a1!=resM and a2!=resN:
        messagebox.showinfo("Error", "wrong size of Matrix to add!")
        return 0
    res = []
    for i in range(resN):
        res.append([])
        for j in range(resM):
            temp = ''
            temp += '((' + str(A[i][j]) + ')-(' + str(B[i][j]) + '))'
            res[i].append(temp)
    return res

def Matplus(A):
    inds = sp.Matrix(A).rref()

    matC = []
    k = 0
    for i in range(len(A)):
        if i in inds[1]:
            matC.append([])
            for j in range(len(A[0])):
                matC[k].append(A[i][j])
            k += 1

    matCplus = np.transpose(matC).dot((np.linalg.inv(np.dot(matC, np.transpose(matC)))))
    matB = np.dot(A, matCplus)
    matBplus = np.linalg.inv(matB.transpose().dot(matB)).dot(matB.transpose())
    matAplus = matCplus.dot(matBplus)
    return np.round(matAplus,3)

def Result():
    global additionalRes
    global matrixA
    global vB
    global T0
    global T1
    global vX
    global str000
    global odn
    AAt = matrixMulti(matrixA,np.transpose(matrixA))
    P1 = []
    i =0
    j = 0
    try:
        for i in range(len(AAt)):
            P1.append([])
            for j in range(len(AAt[0])):
                str000 = AAt[i][j]
                P1[i].append(round(integrate.quad(parseEval,T0,T1)[0],3))
    except:
        messagebox.showinfo("Error!", "error in P1("+str(i)+";"+str(j)+")!")

        return vX
    try:
        P1Plus = Matplus(P1)
    except:
        messagebox.showinfo("Error!", " P1 has no reverse matrix!")
    str000 = ''
    Av = []
    mTemp = matrixMulti(matrixA,vt)
    i = 0
    j = 0
    try:
        for i in range(len(mTemp)):
            Av.append([])
            for j in range(len(mTemp[0])):
                str000 = mTemp[i][j]
                Av[i].append(round(integrate.quad(parseEval,T0,T1)[0],3))
    except:
        messagebox.showinfo("Error!", "error in Av("+str(i)+";"+str(j)+")!")
        return
    odn = np.linalg.det(P1)


    fstM = matrixMulti(np.transpose(matrixA),P1Plus)
    fstM = matrixMulti(fstM,vB)
    sndM = matrixMulti(np.transpose(matrixA),P1Plus)
    sndM = matrixMulti(sndM,Av)

    result = matrixAdd(fstM,vt)
    result = matrixMinus(result,sndM)


    acu = matrixMinus(matrixMulti(np.transpose(vB),vB),matrixMulti(matrixMulti(matrixMulti(np.transpose(vB),P1),P1Plus),vB))

    additionalRes = fstM
    vX = copy.deepcopy(result)
    print(vX)
    resstr =''

    if odn<=0:
        resstr+="the solution isn`t clear\n"
        str000 = acu[0][0]
        for i in range(len(vX)):
            resstr+=vX[i][0]+'\n'
    else:
        resstr += "the solution is clear\n"
        str000 = acu[0][0]
        for i in range(len(vX)):
            resstr += additionalRes[i][0] + '\n'
    messagebox.showinfo("x(t)=",resstr+'\n'+'accuracy = '+str(parseEval(0)))
    str000=''
    return vX

def check_answ():
    global matrixA
    global matrixB
    global T1
    global vX
    global vB
    global str000
    global T0
    global odn
    global additionalRes
    matAX = matrixMulti(matrixA,vX)
    matAX2 = matrixMulti(matrixA, additionalRes)
    vCheck = copy.deepcopy(vB)
    resstr =''
    for i in range(len(vB)):
        str000 = matAX[i][0]
        vCheck = round(integrate.quad(parseEval,T0,T1)[0],3)
        resstr+=str(vCheck)+'\n'
    resstr = ''

    for i in range(len(vB)):
        str000 = matAX2[i][0]
        vCheck = round(integrate.quad(parseEval, T0, T1)[0], 3)
        resstr += str(vCheck) + '\n'
    if odn <=0:
        messagebox.showinfo("Check", resstr + '\n')
    else:
        messagebox.showinfo("Check2",resstr+'\n')



#----------------------------buttons----------------------------
def click_main():
    global m
    global n
    show_message()

def ExitC():
    wind1 = messagebox.askquestion("Bye?!", message="Do you really want to exit?!")
    if wind1 == 'no':
        return
    exit()
def show_message():
    #messagebox.showinfo("m,n", "m = "+str(message.get())+"; n = "+str(message1.get()))
    global m
    global n
    global matrixToDel
    wind1 = messagebox.askquestion("Attention!", message="previous data will be cleared!")
    if wind1 == 'no':
        return
    for i in range(n):
        if vectorToDel[i] != []:
            vectorToDel[i].destroy()
        for j in range(m):
            if matrixToDel[i][j] != []:
                matrixToDel[i][j].destroy()
    m = message.get()
    n = message1.get()
    show_matrix_input()


def click_enter():
    global matrixA
    global matrixToDel
    i = imessage.get()
    j = jmessage.get()
    func = funcmess.get()
    if not parse0(func):
        messagebox.showinfo("Error",func+" is not a function!")
        return
    if i<0 or i >= n:
        messagebox.showinfo("Error",str(i) + " is out of range [0,m) for i! m = "+str(n))
        return
    if j < 0 or j >= m:
        messagebox.showinfo("Error",str(j) + " is out of range [0,m) for j! n = "+str(m))
        return
    if matrixToDel[i][j] != []:
        matrixToDel[i][j].destroy()
    matrixA[i][j]=func
    if(i<10)and j<10:
        tl = Label(text=func,bg='red')
        tl.place(x=j * 100 + 100, width=100, height=20, y=120 + i * 30)
        matrixToDel[i][j]=tl

    return
def click_enter2():
    global vB
    global vectorToDel

    i = imessage2.get()

    bi = Bmess.get()

    if i < 0 or i >= n:
        messagebox.showinfo("Error", str(i) + " is out of range [0,m) for i! m = " + str(n))
        return

    if vectorToDel[i] != []:
        vectorToDel[i].destroy()
    vB[i][0] = bi
    if (i < 10) :
        tl = Label(text=str(bi), bg='red')
        tl.place(x=i * 100 + 100, width=100, height=20, y=120 + n*30+80)
        vectorToDel[i] = tl
    return

    return
def click_clear():
    global clear0
    wind1 = messagebox.askquestion("Attention!", message="Do you really want to clear all?!")
    if wind1 == 'no':
        return

    for i in range(n):
        if vectorToDel[i] != []:
            vectorToDel[i].destroy()
        for j in range(m):
            if matrixToDel[i][j] != []:
                matrixToDel[i][j].destroy()
    for el in clear0:
        el.destroy()
def click_showij():
    global matrixA
    global n
    global m
    i = imessage.get()
    j = jmessage.get()
    if 0<=j<m and 0<=i<n:
        messagebox.showinfo("Show", "a("+str(i)+";"+str(j) + ") = " + str(matrixA[i][j]))
    else:
        messagebox.showinfo("Error", str(i)+" or "+str(j) + " isnt in matrix! n = " + str(m)+"; m = "+str(n))
def click_showi():
    global vB
    global n
    global m
    i = imessage2.get()

    if 0<=i<n:
        messagebox.showinfo("Show", "b("+str(i)+ ") = " + str(vB[i]))
    else:
        messagebox.showinfo("Error", str(i) + " isnt in vector! m = " + str(n))

def click_TT():
    global T0
    global T1
    T0 = T0message.get()
    T1 = T1message.get()

def click_result():
    Result()

def show_matrix_input():
    global clear0
    global matrixA
    global matrixToDel
    global vB
    global vX
    global vt
    global vectorToDel



    matrixA = []
    vB =[]
    vX=[]
    vt=[]
    matrixToDel=[[]]
    for i in range(n):
        matrixToDel.append([])
        matrixA.append([])
        vB.append([0])
        vectorToDel.append([])
        for j in range(m):
            matrixToDel[i].append([])
            matrixA[i].append('0')
    for i in range(m):
        vX.append([])
        vt.append(['pow(t,2)+1'])

    for el in clear0:
        el.destroy()
    btn1 = Button(text="Enter (i;j) element for matrix A(t)", command=click_enter)
    btn1.place(x=0, y=30, width=275, height=26)
    clear0.append(btn1)

    btn2 = Button(text="Enter (i) element for vector b", command=click_enter2)
    btn2.place(x=0, y=60, width=275, height=26)
    clear0.append(btn2)
    iL2 = Label(text="i = ")
    iL2.place(x=300, y=62, width=30, height=20)
    clear0.append(iL2)
    iEnter2 = Entry(textvariable=imessage2)
    iEnter2.place(x=370, y=73, width=75, height=20, anchor="c")
    clear0.append(iEnter2)
    FL2 = Label(text="bi = ")
    FL2.place(x=550, y=62, width=30, height=20)
    clear0.append(FL2)
    BEnter = Entry(textvariable=Bmess)
    BEnter.place(x=670, y=73, width=150, height=20, anchor="c")
    clear0.append(BEnter)
    btnShow2 = Button(text="Show b(i)",command=click_showi)
    btnShow2.place(x=750, y=60, width=70, height=20)
    clear0.append(btnShow2)
    btnShow1 = Button(text="Show a(i;j)", command=click_showij)
    btnShow1.place(x=750, y=30, width=70, height=20)
    clear0.append(btnShow1)
    btnRes = Button(text="Result!", command=click_result)
    btnRes.place(x=850, y=30, width=60, height=20)
    clear0.append(btnRes)
    btnResCheck = Button(text="Check!", command=check_answ)
    btnResCheck.place(x=850, y=60, width=60, height=20)
    clear0.append(btnResCheck)

    iL = Label(text="i = ")
    iL.place(x=300, y=32, width=30, height=20)
    clear0.append(iL)

    iEnter = Entry(textvariable=imessage)
    iEnter.place(x=370, y=43, width=75, height=20, anchor="c")
    clear0.append(iEnter)

    jL = Label(text="j = ")
    jL.place(x=420, y=32, width=30, height=20)
    clear0.append(jL)

    jEnter = Entry(textvariable=jmessage)
    jEnter.place(x=490, y=43, width=75, height=20, anchor="c")
    clear0.append(jEnter)

    FL = Label(text="F(t) = ")
    FL.place(x=550, y=32, width=30, height=20)
    clear0.append(FL)

    FEnter = Entry(textvariable=funcmess)
    FEnter.place(x=670, y=43, width=150, height=20, anchor="c")
    clear0.append(FEnter)
    tl2 = Label(text='A(t)')
    tl2.place(x=30, width=50, height=20, y=100)
    clear0.append(tl2)
    for i in range (n):
        if i == 10:
            break
        tl = Label(text=str(i),bg='red')
        tl.place(x=30, width=50, height=20, y=120 + i * 30)
        clear0.append(tl)
    for j in range(m):
        if j == 10:
            break
        tl = Label(text=str(j),bg='red')
        tl.place(x=j * 100 + 100, width=100, height=20, y=100)
        clear0.append(tl)
    tl3 = Label(text='b')
    tl3.place(x=30, width=50, height=20, y=120+30*n+50)
    clear0.append(tl3)
    for b in range(n):
        if b == 10:
            break
        tl = Label(text=str(b), bg='red')
        tl.place(x=b * 100 + 100, width=100, height=20, y=120+30*n+50)
        clear0.append(tl)

#------------------------main-------------------------------

root = Tk()

additionalRes =[]

root.attributes("-fullscreen",True)

root.title("Лаб2")
root.geometry("700x500")
message = IntVar()
message1 = IntVar()
imessage = IntVar()
imessage2 = IntVar()
jmessage = IntVar()
funcmess = StringVar()
T0message = DoubleVar()
T1message = DoubleVar()
Bmess = DoubleVar()

matrixToDel=[[]]
vectorToDel = []
clear0 = []
n =0
m = 0
T0 = 0
T1 = 0
matrixA = []
vB = []
vX =[]
vt = []
str000 = ''
btn = Button(text="Enter m,n",command=click_main)
btn.place(x=0,y=0,width=75,height=26)
btnEnd = Button(text="Close",foreground='red',command=ExitC)
btnEnd.place(relx=.95,y=0,relwidth=.05)

lm = Label(text="m = ")
lm.place(x=80,y=2,width=30,height=20)

message_entry1 = Entry(textvariable=message1)
message_entry1.place(x=150, y=13,width=75,height=20, anchor="c")

lm1 = Label(text="n = ")
lm1.place(x=180,y=2,width=30,height=20)

message_entry1 = Entry(textvariable=message)
message_entry1.place(x=250, y=13,width=75,height=20, anchor="c")

lm2 = Label(text="T0 = ")
lm2.place(x=300,y=2,width=30,height=20)

message_entry1 = Entry(textvariable=T0message)
message_entry1.place(x=370, y=13,width=75,height=20, anchor="c")

lm13 = Label(text="T1 = ")
lm13.place(x=400,y=2,width=30,height=20)

message_entry1 = Entry(textvariable=T1message)
message_entry1.place(x=470, y=13,width=75,height=20, anchor="c")

btn = Button(text="Enter T0,T1",command=click_TT)
btn.place(x=550,y=0,width=75,height=26)

btn = Button(text="Clear",command=click_clear)
btn.place(relx=.9,y=0,relwidth=.05,height=26)

root.mainloop()


