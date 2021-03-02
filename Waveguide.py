from numpy import (sin,cos,tan,exp,sqrt,pi,arange,meshgrid,linspace)
from tkinter import (Tk,Frame,Button,Entry,Label,StringVar)
from tkinter.ttk import Combobox
from matplotlib import use,rcParams
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.pyplot import (clf,figure,contour,subplot,xlabel,ylabel,subplots_adjust,title)
from matplotlib.ticker import LogLocator
from sys import exit as sys_exit

class WaveG:
    def __init__(self, Moda, n, m, lamb, a, b):
        self.Moda = Moda
        self.n = n
        self.m = m
        self.lamb = lamb
        self.a = a
        self.b = b
        self.h = 0.01 #шаг сетки

    #считаем параметры
    def kappa(self):
        sp = 3 #скорость света в см/сек*10^(10)
        omeg = (2*pi*sp)/self.lamb #частота
        kp = sqrt(((pi*self.n/self.a)**2)+((pi*self.m/self.b)**2)) #каппа
        kpX = pi*self.n/self.a #каппаX
        kpY = pi*self.m/self.b #каппаY
        hh = 2*pi/self.lamb #h
        omegK = sp*kp
        return omeg, kp, kpX, kpY, hh, omegK

    #уравнения задаем
    def makeTEnm(xx,yy,time,omeg,kp,kpX,kpY,hh,self):
        #TE Hxy
        x = arange(self.h, self.a, self.h)
        y = arange(self.h, self.b, self.h)
        xgrid, ygrid = meshgrid(x, y)
        zgrid = ((abs(sin(kpX*xgrid))**((1/kpX)**2))/(abs(sin(kpY*ygrid))**((1/kpY)**2)))*cos(omeg*time)
        #TE Hxz
        x1 = arange(self.h, self.a, self.h)
        z1 = arange(self.h, self.lamb, self.h)
        xgrid1, zgrid1 = meshgrid(x1, z1)
        ygrid1 = abs(sin(kpX*xgrid1))**((kp/kpX)**2)*abs(cos(omeg*time-hh*zgrid1))*cos(kpY*yy)
        #TE Hzy
        z2 = arange(self.h, self.lamb, self.h)
        y2 = arange(self.h, self.b, self.h)
        zgrid2, ygrid2 = meshgrid(z2, y2)
        xgrid2 = abs(sin(kpY*ygrid2))**((kp/kpY)**2)*abs(cos(omeg*time-hh*zgrid2))*cos(kpX*xx)
        #TE Exy
        x3 = arange(self.h, self.a, self.h)
        y3 = arange(self.h, self.b, self.h)
        xgrid3, ygrid3 = meshgrid(x3, y3)
        zgrid3 = abs(cos(kpY*ygrid3))*abs(cos(kpX*xgrid3))*cos(omeg*time)
        #TE Exz
        x4 = arange(self.h, self.a, self.h)
        z4 = arange(self.h, self.lamb, self.h)
        xgrid4, zgrid4 = meshgrid(x4, z4)
        ygrid4 = (abs(cos(kpX*xgrid4))/(exp(kpY*tan(kpY*yy)*zgrid4)))*cos(omeg*time-hh*zgrid4)
        #TE Ezy
        z5 = arange(self.h, self.lamb, self.h)
        y5 = arange(self.h, self.b, self.h)
        zgrid5, ygrid5 = meshgrid(z5, y5)
        xgrid5 = (abs(cos(kpY*ygrid5))/(exp(kpX*tan(kpX*xx)*zgrid5)))*cos(omeg*time-hh*zgrid5)
        return xgrid, ygrid, zgrid, xgrid1, ygrid1, zgrid1, xgrid2, ygrid2, zgrid2, xgrid3, ygrid3, zgrid3, xgrid4, ygrid4, zgrid4, xgrid5, ygrid5, zgrid5

    def makeTMnm(xx,yy,time,omeg,kp,kpX,kpY,hh,self):
        #TM Hxy
        x = arange(self.h, self.a, self.h)
        y = arange(self.h, self.b, self.h)
        xgrid, ygrid = meshgrid(x, y)
        zgrid = abs(sin(kpX*xgrid))*abs(sin(kpY*ygrid))*cos(omeg*time)
        #TM Hxz
        x1 = arange(self.h, self.a, self.h)
        z1 = arange(self.h, self.lamb, self.h)
        xgrid1, zgrid1 = meshgrid(x1, z1)
        ygrid1 = abs(sin(kpX*xgrid1))*exp((1/tan(kpY*yy))*kpY*zgrid1)*cos(omeg*time-hh*zgrid1)
        #TM Hzy
        z2 = arange(self.h, self.lamb, self.h)
        y2 = arange(self.h, self.b, self.h)
        zgrid2, ygrid2 = meshgrid(z2, y2)
        xgrid2 = abs(sin(kpY*ygrid2))*exp((1/tan(kpX*xx))*kpX*zgrid2)*cos(omeg*time-hh*zgrid2)
        #TM Exy
        x3 = arange(self.h, self.a, self.h)
        y3 = arange(self.h, self.b, self.h)
        xgrid3, ygrid3 = meshgrid(x3, y3)
        zgrid3 = ((abs(cos(kpY*ygrid3))**(1/(kpY**2)))/(abs(cos(kpX*xgrid3))**(1/(kpX**2))))*cos(omeg*time)
        #TM Exz
        x4 = arange(self.h, self.a, self.h)
        z4 = arange(self.h, self.lamb, self.h)
        xgrid4, zgrid4 = meshgrid(x4, z4)
        ygrid4 = (abs(cos(kpX*xgrid4))**((kp/kpX)**2))*abs(cos(omeg*time-hh*zgrid4))*sin(kpY*yy)
        #TM Ezy
        z5 = arange(self.h, self.lamb, self.h)
        y5 = arange(self.h, self.b, self.h)
        zgrid5, ygrid5 = meshgrid(z5, y5)
        xgrid5 = (abs(cos(kpY*ygrid5))**((kp/kpY)**2))*abs(cos(omeg*time-hh*zgrid5))*sin(kpX*xx)
        return xgrid, ygrid, zgrid, xgrid1, ygrid1, zgrid1, xgrid2, ygrid2, zgrid2, xgrid3, ygrid3, zgrid3, xgrid4, ygrid4, zgrid4, xgrid5, ygrid5, zgrid5

    def makeTEn0(xx,yy,time,omeg,kp,kpX,kpY,hh,self):
        #TE Hxy
        x = arange(self.h, self.a, self.h)
        y = arange(self.h, self.b, self.h)
        xgrid, ygrid = meshgrid(x, y)
        zgrid = abs(sin((pi/self.b)*ygrid)*cos(omeg*time))
        #TE Hxz
        x1 = arange(self.h, self.a, self.h)
        z1 = arange(self.h, self.lamb, self.h)
        xgrid1, zgrid1 = meshgrid(x1, z1)
        ygrid1 = abs(sin(kpX*xgrid1))*abs(cos(omeg*time-hh*zgrid1))
        #TE Hzy
        z2 = arange(self.h, self.lamb, self.h)
        y2 = arange(self.h, self.b, self.h)
        zgrid2, ygrid2 = meshgrid(z2, y2)
        xgrid2 = abs(cos(omeg*time-hh*zgrid2))*exp(kpX*ygrid2/tan(kpX*xx))
        #TE Exy
        x3 = arange(self.h, self.a, self.h)
        y3 = arange(self.h, self.b, self.h)
        xgrid3, ygrid3 = meshgrid(x3, y3)
        zgrid3 = cos(kpX*xgrid3)*sin(omeg*time)
        #TE Exz
        x4 = arange(self.h, self.a, self.h)
        z4 = arange(self.h, self.lamb, self.h)
        xgrid4, zgrid4 = meshgrid(x4, z4)
        ygrid4 = sin(kpX*xgrid4)*sin(omeg*time-hh*zgrid4)
        #TE Ezy
        z5 = arange(self.h, self.lamb, self.h)
        y5 = arange(self.h, self.b, self.h)
        zgrid5, ygrid5 = meshgrid(z5, y5)
        xgrid5 = sin(kpX*xx)*sin(omeg*time-hh*zgrid5)
        return xgrid, ygrid, zgrid, xgrid1, ygrid1, zgrid1, xgrid2, ygrid2, zgrid2, xgrid3, ygrid3, zgrid3, xgrid4, ygrid4, zgrid4, xgrid5, ygrid5, zgrid5

    def makeTE0m(xx,yy,time,omeg,kp,kpX,kpY,hh,self):
        #TE Hxy
        x = arange(self.h, self.a, self.h)
        y = arange(self.h, self.b, self.h)
        xgrid, ygrid = meshgrid(x, y)
        zgrid = abs(sin((pi/self.a)*xgrid)*cos(omeg*time))
        #TE Hxz
        x1 = arange(self.h, self.a, self.h)
        z1 = arange(self.h, self.lamb, self.h)
        xgrid1, zgrid1 = meshgrid(x1, z1)
        ygrid1 = abs(cos(omeg*time-hh*zgrid1))*exp(kpY*xgrid1/tan(kpY*yy))
        #TE Hzy
        z2 = arange(self.h, self.lamb, self.h)
        y2 = arange(self.h, self.b, self.h)
        zgrid2, ygrid2 = meshgrid(z2, y2)
        xgrid2 = abs(sin(kpY*ygrid2))*abs(cos(omeg*time-hh*zgrid2))
        #TE Exy
        x3 = arange(self.h, self.a, self.h)
        y3 = arange(self.h, self.b, self.h)
        xgrid3, ygrid3 = meshgrid(x3, y3)
        zgrid3 = cos(kpY*ygrid3)*sin(omeg*time)
        #TE Exz
        x4 = arange(self.h, self.a, self.h)
        z4 = arange(self.h, self.lamb, self.h)
        xgrid4, zgrid4 = meshgrid(x4, z4)
        ygrid4 = sin(kpY*yy)*sin(omeg*time-hh*zgrid4)
        #TE Ezy
        z5 = arange(self.h, self.lamb, self.h)
        y5 = arange(self.h, self.b, self.h)
        zgrid5, ygrid5 = meshgrid(z5, y5)
        xgrid5 = sin(kpY*ygrid5)*sin(omeg*time-hh*zgrid5)
        return xgrid, ygrid, zgrid, xgrid1, ygrid1, zgrid1, xgrid2, ygrid2, zgrid2, xgrid3, ygrid3, zgrid3, xgrid4, ygrid4, zgrid4, xgrid5, ygrid5, zgrid5

    #проверка моды и строитель графиков
    def stroit(self,time,omeg,kp,kpX,kpY,hh,omegK):
        rcParams['contour.negative_linestyle'] = 'solid'
        #Составление массива линий уровня !Можно сделать динамическим
        def arr2(arrayAdd=[]):
            step=0
            arrayValues=[]
            plus=0
            minus=0
            while step<50:
                plus+=1/(2**step)
                minus-=1/(2**step)
                arrayValues.append(plus-1)
                arrayValues.append(minus+1)
                step+=1
            arrayValues.append(0.15)
            arrayValues.append(-0.15)
            arrayValues1=sorted(arrayValues)
            del arrayValues1[51]
            del arrayValues1[50]
            arrayValues1+=arrayAdd
            arrayValues2=sorted(arrayValues1)
            return arrayValues2

        def arrS(arrayAdd=[]):
            stepX=pi/10
            XX=0
            arrayValues=[]
            while XX<=pi/2:
                arrayValues.append(sin(XX))
                XX+=stepX
            arrayValues+=arrayAdd
            return arrayValues

        #метки для регистрации ошибок
        N1=1
        N2=1
        M1=0
        if omeg>omegK:
            if self.Moda==0:
                if self.n==0 and self.m==0:
                    M1=1
                    return M1
                elif self.n==0 and self.m!=0:
                    if self.m%2==0 and self.m!=0:
                        N2=2 #четное ли б, если да то делим б на 2 для выбора сечения волновода в нужном месте
                    x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5 = WaveG.makeTE0m(self.a,self.b/(2*self.m),time,omeg,kp,kpX,kpY,hh,self)
                    subplot(2,4,1)
                    C1=contour(x, y, z, arrS(), linewidths=1.8, colors='yellow')
                    xlabel("x")
                    ylabel("y")
                    title("Hxy")
                    subplot(2,4,5)
                    C2=contour(x1, z1, y1, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("z")
                    title("Hxz")
                    subplot(2,4,2)
                    C3=contour(z2, y2, x2, arr2(), linewidths=1.8)
                    xlabel("z")
                    ylabel("y")
                    title("Hzy")
                    subplot(2,4,3)
                    C4=contour(x3, y3, z3,linspace(-1,1,20), linewidths=1.8, colors='yellow')
                    xlabel("x")
                    ylabel("y")
                    title("Exy")
                    subplot(2,4,7)
                    C5=contour(x4, z4, y4, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("z")
                    title("Exz")
                    subplot(2,4,4)
                    C6=contour(z5, y5, x5, arr2(), linewidths=1.8)
                    xlabel("z")
                    ylabel("y")
                    title("Ezy")
                    subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.99, hspace=0.25, wspace=0.25)
                elif self.n!=0 and self.m==0:
                    if self.n%2==0 and self.n!=0:
                        N1=2 #четное ли n, если да то делим а на 2 для выбора сечения волновода в нужном месте
                    x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5 = WaveG.makeTEn0(self.a/(2*self.n),self.b,time,omeg,kp,kpX,kpY,hh,self)
                    subplot(2,4,1)
                    C1=contour(x, y, z, arrS(), linewidths=1.8, colors='yellow')
                    xlabel("x")
                    ylabel("y")
                    title("Hxy")
                    subplot(2,4,5)
                    C2=contour(x1, z1, y1, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("z")
                    title("Hxz")
                    subplot(2,4,2)
                    C3=contour(z2, y2, x2, arr2(), linewidths=1.8)
                    xlabel("z")
                    ylabel("y")
                    title("Hzy")
                    subplot(2,4,3)
                    C4=contour(x3, y3, z3,linspace(-1,1,20), linewidths=1.8, colors='yellow')
                    xlabel("x")
                    ylabel("y")
                    title("Exy")
                    subplot(2,4,7)
                    C5=contour(x4, z4, y4, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("z")
                    title("Exz")
                    subplot(2,4,4)
                    C6=contour(z5, y5, x5, arr2(), linewidths=1.8)
                    xlabel("z")
                    ylabel("y")
                    title("Ezy")
                    subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.99, hspace=0.25, wspace=0.25)
                else:
                    x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5 = WaveG.makeTEnm(self.a,self.b,time,omeg,kp,kpX,kpY,hh,self)
                    subplot(2,4,1)
                    C1=contour(x, y, z,  linspace(-2,2,40), linewidths=1.8)
                    xlabel("x")
                    ylabel("y")
                    title("Hxy")
                    subplot(2,4,5)
                    C2=contour(x1, z1, y1, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("z")
                    title("Hxz")
                    subplot(2,4,2)
                    C3=contour(z2, y2, x2, arr2(), linewidths=1.8)
                    xlabel("z")
                    ylabel("y")
                    title("Hzy")
                    subplot(2,4,3)
                    C4=contour(x3, y3, z3, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("y")
                    title("Exy")
                    subplot(2,4,7)
                    C5=contour(x4, z4, y4, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("z")
                    title("Exz")
                    subplot(2,4,4)
                    C6=contour(z5, y5, x5, arr2(), linewidths=1.8)
                    xlabel("z")
                    ylabel("y")
                    title("Ezy")
                    subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.99, hspace=0.25, wspace=0.25)
            elif self.Moda==1:
                if self.n!=0 and self.m!=0:
                    if self.n%2==0 and self.n!=0:
                        N1=2 #четное ли n, если да то делим а на 2 для выбора сечения волновода в нужном месте
                    if self.m%2==0 and self.m!=0:
                        N2=2 #четное ли б, если да то делим б на 2 для выбора сечения волновода в нужном месте
                    x, y, z, x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4, x5, y5, z5 = WaveG.makeTMnm(self.a/(2*self.n),self.b/(2*self.m),time,omeg,kp,kpX,kpY,hh,self)
                    subplot(2,4,1)
                    C1=contour(x, y, z, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("y")
                    title("Hxy")
                    subplot(2,4,5)
                    C2=contour(x1, z1, y1, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("z")
                    title("Hxz")
                    subplot(2,4,2)
                    C3=contour(z2, y2, x2, arr2(), linewidths=1.8)
                    xlabel("z")
                    ylabel("y")
                    title("Hzy")
                    subplot(2,4,3)
                    if abs(cos(omeg*time))>=0.0001:
                        C4=contour(x3, y3, z3, linspace(-1,1,40), linewidths=1.8)
                    else:
                        NN=1
                    xlabel("x")
                    ylabel("y")
                    title("Exy")
                    subplot(2,4,7)
                    C5=contour(x4, z4, y4, arr2(), linewidths=1.8)
                    xlabel("x")
                    ylabel("z")
                    title("Exz")
                    subplot(2,4,4)
                    C6=contour(z5, y5, x5, arr2(), linewidths=1.8)
                    xlabel("z")
                    ylabel("y")
                    title("Ezy")
                    subplots_adjust(top=0.9, bottom=0.1, left=0.1, right=0.99, hspace=0.25, wspace=0.25)
                else:
                    M1=1
                    return M1
        else:
            M1=2
            if self.Moda==0:
                if self.n==0 and self.m==0:
                    M1+=1
            elif self.Moda==1:
                if self.n==0 or self.m==0:
                    M1+=1
            return M1

        #строитель графиков
        #[i*sin(omeg*time) for i in arr2(arrAddEE)] умножениеэлементов массива на число
        #abs(linspace(-5,5,20)) толщина линий регулируется
        #extend='both' область вне диапозона закрашивает
        #linspace(-1,1,20)
        return M1

class Main(Frame):
    def __init__(self, root):
        super().__init__(root)

#проверка что модуль запущен самостоятельно пользователем
if __name__ == "__main__":
    #Функция вызывающая график с параметрами задаными в главном окне (ВНулевойМоментВремениПоУмолчанию #время в сек*10^(-10))
    def Vvod(time):
        label8=Label(frame2,text=' ',width=30,height=2)
        #проверка на то число ли n и m?
        if ((message1.get()).isdigit() and (message2.get()).isdigit() and (message3.get()).replace('.','',1).isdigit() and (message4.get()).replace('.','',1).isdigit() and\
        (message5.get()).replace('.','',1).isdigit() and float(message3.get())!=0 and float(message4.get())!=0 and float(message5.get())!=0):
            Moda = combobox1.current()
            plot1 = WaveG(Moda,float(message1.get()),float(message2.get()),float(message3.get()),float(message4.get()),float(message5.get()))
            omeg,kp,kpX,kpY,hh,omegK = plot1.kappa() #получаем данные
            clf() #ОчиститьТекущийРисунок
            M1 = plot1.stroit(time,omeg,kp,kpX,kpY,hh,omegK) #ВызовФункцииПостроения
            fig.canvas.draw() #РисуемПоверхСтарого
            label7=Label(frame2,text='f={0}ГГц, fкр={1}ГГц'.format(round(omeg*10/(2*pi),2),round(omegK*10/(2*pi),2)),width=30,height=1) #Вывод частот в окно с округлением до сотых
            label7.grid(row=0,column=1)
            if M1==0:
                label8.config(text='Мода есть')
            elif M1==1:
                label8.config(text='Моды нет')
            elif M1==2:
                label8.config(text='Частота меньше критической')
            else:
                label8.config(text='Моды нет \nЧастота меньше критической')
            #еше заголовок,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
        else:
            label8.config(text='Неверно заполнено')
        label8.grid(row=1,column=1)
    #время вперед
    def TimeS():
        Moda = combobox1.current()
        plot1 = WaveG(Moda,float(message1.get()),float(message2.get()),float(message3.get()),float(message4.get()),float(message5.get()))
        omeg,kp,kpX,kpY,hh,omegK = plot1.kappa()
        global time
        time+=0.1*pi/omeg
        if time<=0.001:
            time=0
    #сброс времени
    def TimeRest():
        global time
        time=0

    #Создание главного окна
    use('TkAgg') #интерактивный бэкэнд
    root = Tk()
    root.title("Waveguide")
    root.geometry("1700x800") #размер окна программы
    root.resizable(False,False) #запретить изменение размера
    app = Main(root)
    app.pack()

    #Добавление графика в окно
    fig = figure(figsize=(16,8), dpi=90) #dpi регулирует размер окна с графиком, figsize соотношение сторон окна с графиком
    canvas = FigureCanvasTkAgg(fig, master=root)
    plot_widget = canvas.get_tk_widget()
    plot_widget.pack(side='right')

    #Создание кнопок и полей ввода
    frame1=Frame()
    frame2=Frame()

    label1=Label(frame1,text=u'Мода')
    combobox1=Combobox(frame1, values = ["TE","TM"], state='readonly')
    combobox1.current(0)  #вариант по умолчанию
    label1.grid(row=0,column=0)
    combobox1.grid(row=0,column=1)

    label2=Label(frame1,text=u'n',width=6)
    label2.grid(row=1,column=0)
    label3=Label(frame1,text=u'm',width=6)
    label3.grid(row=1,column=1)

    message1=StringVar()
    entry1=Entry(frame1,textvariable=message1,width=6)
    entry1.grid(row=2,column=0)
    message2=StringVar()
    entry2=Entry(frame1,textvariable=message2,width=6)
    entry2.grid(row=2,column=1)

    label4=Label(frame1,text=u'lam (см):',width=6)
    message3=StringVar()
    entry3=Entry(frame1,textvariable=message3,width=6)
    label4.grid(row=3,column=0)
    entry3.grid(row=3,column=1)

    label5=Label(frame1,text=u'a (см):',width=6)
    message4=StringVar()
    entry4=Entry(frame1,textvariable=message4,width=6)
    label5.grid(row=4,column=0)
    entry4.grid(row=4,column=1)

    label6=Label(frame1,text=u'b (см):',width=6)
    message5=StringVar()
    entry5=Entry(frame1,textvariable=message5,width=6)
    label6.grid(row=5,column=0)
    entry5.grid(row=5,column=1)

    time = 0
    button2=Button(root,text=u'"Шаг времени вперед"',command=lambda:[TimeS(),Vvod(time)]) #лямбда-выражение
    button2.pack(side='bottom')

    button1=Button(root,text=u'Посторить график',command=lambda:[TimeRest(),Vvod(time)])
    button1.pack(side='bottom')

    frame1.place(x=10, y=20)
    frame2.place(x=10, y=200)

    root.protocol("WM_DELETE_WINDOW", sys_exit)

    root.mainloop()
