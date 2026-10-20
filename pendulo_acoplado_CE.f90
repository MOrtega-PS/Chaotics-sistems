MODULE constantes 
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::ki=SELECTED_REAL_KIND(15,300)
END MODULE constantes

Module varglob
  use constantes
  IMPLICIT NONE
  real(ki):: ke, q1, q2, l1,l2,x01,x02,g,m1,m2
end module varglob

program pendulos_cargados
   use constantes
  use varglob
  implicit none
  integer :: i,n,auxp,auxc
  real(ki), allocatable:: os1(:),ws1(:),os2(:),ws2(:)
  real(ki),allocatable:: xs1(:),ys1(:),xs2(:),ys2(:),p1(:),p2(:)
  real(ki):: a,b,h,ts,aux,En, k1(4),k2(4),k3(4),k4(4), veact(4),vesig(4)
  real(ki)::tini,tfin,t_cal

  !!!!diccionario de variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !a,b limites del intervalo
  !h tamaño de paso
  !n numero de subintervalos
  !i variable de conteo
  !os1 os2 Puntos de la solucion para el angulo
  !ws1 ws2 Puntos de la solucion para la velocidad angular
  !k Valores de la k para runge kuta
  !aux variable auxiliar en caso de necesitarla
  !xs Posicion en el eje x
  !ys Posicion en el eje y
  !En energia del sistema
  !ke constante de Coulomb (1/(4 pi eps0))
  !q1,q2 cargas de las lentejas
  !g gravedad
  !x01 posicion del primer pivote (-x01,0)
  !x02 posicion del segundo pivote (x02,0)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call cpu_time(tini)
  !!!!!!!!!!!!se le pide al usuario los val iniciales!!!!!!!!!!!!!!!
  a = 0
  b = 10
  n = 10000000
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!! se inicializan variables y arrays!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  h = (b-a)/n!se calcula el paso

  allocate(ws1(n+1)) !w(t) velocidad angular
  allocate(os1(n+1)) !0(t) angulo
  allocate(xs1(n+1))  !posicion en x
  allocate(ys1(n+1))  !posicion en y
  allocate(ws2(n+1)) !w(t) velocidad angular
  allocate(os2(n+1)) !0(t) angulo
  allocate(xs2(n+1))  !posicion en x
  allocate(ys2(n+1))  !posicion en y
  
  !!!!!valores automaticos para pruebas rapidas

  !pendulo 1
  ws1(1) = 0.0_ki
  os1(1) = -0.1_ki
  m1 = 1.5_ki
  x01 = 0.5_ki
  l1 = 1.0_ki

  !pendulo 2
  ws2(1) = 0.0_ki
  os2(1) = 1.0_ki
  m2 = 2.0_ki
  x02 = 0.5_ki
  l2 = 0.9_ki

  !cargas (mismo signo = repulsion, signo distinto = atraccion)
  ke = 8.9875517923e9_ki
  q1 = -2.5e-5_ki
  q2 = 2.0e-5_ki
  g = 9.81_ki 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !se abre el archivo donde se guarda la solucion
  open(unit=2,file ="data_polares_pc",status="replace")
  open(unit=3,file ="data_cartesianas_pc",status="replace")

  !!!!!!!!!!!!!!!!!!!!metodo runge kuta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !incializamos el vector de estado
  veact= [os1(1),os2(1),ws1(1),ws2(1)]
  do i =1,n
     ts=a+(i-1)*h !guardas la t que se "usa"
     !!!!!!calculo de k, vectorial!!!!!!
     k1 = h*dv(ts, veact)
     k2 = h*dv(ts + 0.5_ki*h, veact + 0.5_ki*k1)
     k3 = h*dv(ts + 0.5_ki*h, veact + 0.5_ki*k2)
     k4 = h*dv(ts + h, veact + k3)

     vesig= veact + (k1 + 2.0_ki*k2 + 2.0_ki*k3 + k4) / 6.0_ki
     !se guardan los resultados en los vectores
     os1(i+1) = vesig(1)
     os2(i+1) = vesig(2)
     ws1(i+1) = vesig(3)
     ws2(i+1) = vesig(4)

     veact = vesig

  end do
!!!!!!!!!!!!!!!!!!!!!!!!!fin de rungekuta!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!! Conversion de polares a cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i =1,n+1
     xs1(i) = -x01 + l1*sin(os1(i))
     ys1(i) = -l1*cos(os1(i))
     xs2(i) = x02 + l2*sin(os2(i))
     ys2(i) = -l2*cos(os2(i))
  end do
  
  !!!!!!!!!!! Fin de conversion a cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!escritura de las tabla cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ts = a
  auxc = 5000
  do i =0,auxc
     
     write(3,*) ts, xs1((n/auxc)*i+1),ys1((n/auxc)*i+1),xs2((n/auxc)*i+1),ys2((n/auxc)*i+1)
      ts = ts+h*(n/auxc)
  enddo
!!!!!!!!!!!!!fin de la escritura de las tablas cartesianas!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!escritura de las tablas polares!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ts = a
  auxp= 5000
  do i =0,auxp
     
     write(2,*) ts, os1((n/auxp)*i+1),ws1((n/auxp)*i+1),os2((n/auxp)*i+1),ws2((n/auxp)*i+1)
      ts = ts+h*(n/auxp)
  enddo
!!!!!!!!!!!!!fin de la escritura de las tablas polares!!!!!!!!!!!!!!!!!

  !se cierran los archivos de datos para que gnuplot los lea completos
  close(2)
  close(3)

  !!!!!!!!!!!!!!!!!!!!!!!!!graficacion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!grafica 1!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit= 10, file = "datos.plt", status = "unknown")
  write(10,*)'m = "./data_polares_pc"'!indicaciones del archivo de datos
  write(10,*) 'set title "Posicion angular 1 contra tiempo"'
  write(10,*) 'set autoscale'
  write(10,*) 'unset key'
  write(10,*) 'set xlabel "Tiempo"'
  write(10,*) 'set ylabel "Posicion angular"'
  write(10,*)'plot m using 1:2 with lines'
  close(10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!graficacion!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!grafica 2!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit= 11, file = "datos2.plt", status = "unknown")
  write(11,*)'m = "./data_polares_pc"'!indicaciones del archivo de datos
  write(11,*) 'set title "Posicion agular 2 contra tiempo"'
  write(11,*) 'set autoscale'
  write(11,*) 'unset key'
  write(11,*) 'set xlabel "Tiempo"'
  write(11,*) 'set ylabel "Posicion angular"'
  write(11,*)'plot m using 1:4 with lines'
  close(11)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!grafica 3: espacio fase pendulo 1!!!!!!!!!!!!!!!!!!
  open(unit= 12, file = "datos3.plt", status = "unknown")
  write(12,*)'m = "./data_polares_pc"'
  write(12,*) 'set title "Espacio fase"'
  write(12,*) 'set autoscale'
  write(12,*) 'unset key'
  write(12,*) 'set xlabel "Posicion angular (rad)"'
  write(12,*) 'set ylabel "Velocidad angular (rad/s)"'
  write(12,*)'plot m using 2:3 with lines title "Pendulo 1", m using 4:5 &
     &with lines title "Pendulo 2"'
  close(12)



  call system("gnuplot -p datos.plt")!llamar al sistema para usar gnuplot
  call system("gnuplot -p datos2.plt")
  call system("gnuplot -p datos3.plt")
  call system("gnuplot -p datos4.plt")
  call cpu_time(tfin)
  t_cal = tfin-tini
  print '(A, F10.4, A)', " Tiempo de ejecución del CPU: ", t_cal, " se&
       &gundos."

  
contains

  function dv(t, ve) result(dve)
    use constantes
    use varglob
    implicit none
    real(ki), intent(in) :: t
    real(ki), intent(in) :: ve(4)
    real(ki) :: dve(4)
    real(ki) :: o1, w1,o2,w2,a,b,c,d,d3,kq
    o1 = ve(1)
    o2 = ve(2)
    w1 = ve(3)
    w2 = ve(4)

    kq = ke*q1*q2
    a = o1-o2; b=x01+x02;c=b**2 + l1**2 + l2**2
    d  = sqrt(c-2*l1*l2*cos(a)+2*b*(l2*sin(o2)-l1*sin(o1)))
    d3 = d**3
    dve(1) = w1
    dve(2) = w2
    dve(3) = (-g/l1)*sin(o1) - kq/(m1*l1*d3)*( b*cos(o1) - l2*sin(a) )
    dve(4) = (-g/l2)*sin(o2) + kq/(m2*l2*d3)*( b*cos(o2) - l1*sin(a) )

  end function dv
end program pendulos_cargados
