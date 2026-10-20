MODULE constantes 
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::ki=SELECTED_REAL_KIND(15,300)
END MODULE constantes

Module varglob
  use constantes
  IMPLICIT NONE
  real(ki)::g,m1,r01,r02,kr1,kr2,m2
end module varglob

program pendulo_doble_elastico_v2
  use constantes
  use varglob
  implicit none
  integer :: i,n,aux2,aux3
  real(ki), allocatable:: rs11(:), rs21(:),os11(:),os21(:)
  real(ki), allocatable:: rs12(:), rs22(:),os12(:),os22(:)
  real(ki),allocatable:: xs(:),ys(:),xs2(:),ys2(:)
  real(ki):: a,b,h,ts,aux,En, a1,a2,a3,a4,t1,t2,vt
  real(ki):: k1(8),k2(8),k3(8),k4(8),veact(8),vesig(8)
  real(ki):: r1,r2,w1,w2,o1,o2,vr1,vr2,tini,tfin,t_cal

  !!!!diccionario de variables1!!!!!!!!!!1!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !a,b limites del intervalo
  !h tamaño de paso
  !n numero de subintervalos
  !i variable de conteo
  !rs1 Puntos de la solucion para la velocidad radial
  !rs2 Puntos de la solucion para la posicion radial
  !os1 Puntos de la solucion para el angulo
  !os2 Puntos de la solucion para la velocidad angular
  !k Valores de la k para runge kuta
  !aux, a1, a2, a3, a4, aux2  variables auxiliares
  !xs Posicion en el eje x
  !ys Posicion en el eje y
  !En energia del sistema
  !ve vector de estado del sistema (siguiente y actual)
  !m,b,x matrices auxiliares para resolver el SEL
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!se le pide al usuario los val iniciales!!!!!!!!!!!!!!!
  !print*,"Establezca la a y despues la b"
  !read*,a,b
  !print*,"Establezca el numero de subintervalos"
  call cpu_time(tini)
  !read*,n
  a = 0
  b = 10
  n = 100000000
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!! se inicializan variables y arrays!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  h = (b-a)/n!se calcula el paso

  !pendulo 1
  allocate(rs11(n+1)) !vr(t) velocidad radial
  allocate(rs21(n+1)) !r(t) posicion radial
  allocate(os11(n+1)) !w(t) velocidad angular
  allocate(os21(n+1)) !0(t) angulo
  allocate(xs(n+1))  !posicion en x
  allocate(ys(n+1))  !posicion en y

  !pendulo 2
  allocate(rs12(n+1)) !vr(t) velocidad radial
  allocate(rs22(n+1)) !r(t) posicion radial
  allocate(os12(n+1)) !w(t) velocidad angular
  allocate(os22(n+1)) !0(t) angulo
  allocate(xs2(n+1))  !posicion en x
  allocate(ys2(n+1))  !posicion en y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!se preguntan los val iniciales de las Ecuaciones!!!!!!!!!!!!!!!!!!!!!!
  !print*,"Escriba la informacion del primer pendulo"
  !print*, "Establezca la posicion radial y angular inicial (en ese orden)"
  !read*,aux
  !rs21(1) = aux
  !r01 = aux !posicion angular inicial
  !read*,aux
  !os21(1) = aux
  !print*, "Establezca la velocidad radial y angular inicial (en ese orden)"
  !read*,aux
  !rs11(1) = aux
  !read*,aux
  !os11(1) = aux
  !print*, "Establezca la masa y la constante k del resorte (en ese orden)"
  !read*,m1
  !read*,kr1

  !print*,"Escriba la informacion del segundo pendulo"
  !print*, "Establezca la posicion radial y angular inicial (en ese orden)"
  !read*,aux
  !rs22(1) = aux
  !r02 = aux !posicion angular inicial
  !read*,aux
  !os22(1) = aux
  !print*, "Establezca la velocidad radial y angular inicial (en ese orden)"
  !read*,aux
  !rs12(1) = aux
  !read*,aux
  !os12(1) = aux
  !print*, "Establezca la masa y la constante k del resorte (en ese orden)"
  !read*,m2
  !read*,kr2
  g = 9.81_ki !se establece la gravedad
  !!!!!valores automaticos para pruebas rapidas
  rs21(1) = 1.5_ki
  os21(1) = 0.2_ki
  rs11(1) = 0.0_ki
  os11(1) = 0.0_ki
  r01 = 1.0_ki
  m1 = 1.0_ki
  kr1 = 1000.0_ki

  rs22(1) = 1.0_ki
  os22(1) = 0.1_ki
  rs12(1) = 0.0_ki
  os12(1) = 0.0_ki
  r02 = 1.5_ki
  m2 = 10.0_ki
  kr2 = 1000.0_ki

  ! Asignar valores iniciales al vector de estado veact
  veact = [rs21(1), os21(1), rs22(1), os22(1), rs11(1), os11(1), rs12(1), os12(1)]

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !se abre el archivo donde se guarda la solucion
  open(unit=2,file ="data_polares",status="replace")
  open(unit=3,file ="data_cartesianas",status="replace")
  open(unit=1,file ="debug",status="replace")


 !!!!!!!!!!!!!!!!!!!!metodo runge kuta !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i =1,n+1
     ts=a+(i-1)*h !guardas la t que se "usa"
     !!!!!!calculo de k, vectorial!!!!!!
     k1 = h*dv(ts, veact)
     k2 = h*dv(ts + 0.5_ki*h, veact + 0.5_ki*k1)
     k3 = h*dv(ts + 0.5_ki*h, veact + 0.5_ki*k2)
     k4 = h*dv(ts + h, veact + k3)

     vesig= veact + (k1 + 2.0_ki*k2 + 2.0_ki*k3 + k4) / 6.0_ki

     !se guardan los resultados en los vectores
     rs21(i+1) = vesig(1)
     os21(i+1) = vesig(2)
     rs22(i+1) = vesig(3)
     os22(i+1) = vesig(4)
     rs11(i+1) = vesig(5)
     os11(i+1) = vesig(6)
     rs12(i+1) = vesig(7)
     os12(i+1) = vesig(8)
     
     veact = vesig

  end do
!!!!!!!!!!!!!!!!!!!!!!!!!fin de rungekuta!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!! Conversion de polares a cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do i =1,n+1
     xs(i) = rs21(i)*sin(os21(i))
     ys(i) = -rs21(i)*cos(os21(i))
     xs2(i) = xs(i) + rs22(i)*sin(os22(i))
     ys2(i) = ys(i) -rs22(i)*cos(os22(i))
  end do
  
  !!!!!!!!!!! Fin de conversion a cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!escritura de las tabla cartesianas!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ts = a
  do i =0,1000
     
     write(3,*) ts, xs((n/1000)*i+1),ys((n/1000)*i+1), xs2((n/1000)*i+1),ys2((n/1000)*i+1)
      ts = ts+h*(n/1000)
  enddo
!!!!!!!!!!!!!fin de la escritura de las tablas cartesianas!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!escritura de las tabla polares!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ts = a
  aux2=0
  aux3=1000
  do i =1,aux3
      r1 = rs21(i); o1 = os21(i); r2 = rs22(i); o2 = os22(i)
      vr1 = rs11(i); w1 = os11(i); vr2 = rs12(i); w2 = os12(i)
      a1 = o1-o2
      t1= (1.0_ki/2.0_ki)*m1*(vr1**2 + (r1**2)*(w1**2))
      t2= (1.0_ki/2.0_ki)*m2*(vr2**2 + (r2**2)*(w2**2)+vr1**2+(r1**2)*(w1**2)+2*vr1*vr2*cos(a1) &
           &-2*r1*vr2*w1*sin(a1)+2*r2*vr1*w2*sin(a1)+2*r1*r2*w1*w2*cos(a1))
      vt= m1*g*ys(i)+m2*g*ys2(i)+(1.0_ki/2.0_ki)*kr1*(r1-r01)**2 + (1.0_ki/2.0_ki)*kr2*(r2-r02)**2
      En= t1+t2+vt
      write(2,*) ts, rs21((n/aux3)*i+1),rs22((n/aux3)*i+1), os21((n/aux3)*i+1),os22((n/aux3)*i+1), &
                 & rs11((n/aux3)*i+1),rs12((n/aux3)*i+1), os11((n/aux3)*i+1),os12((n/aux3)*i+1), En
      ts = ts+h*(n/aux3)
     
      !write(2,*) ts, rs21(i+aux2),rs22(i+aux2), os21(i+aux2),os22(i+aux2), &
      !   & rs11(i+aux2),rs12(i+aux2), os11(i+aux2),os12(i+aux2)
      !ts = ts+h
     
  enddo
  !!!!!!!!!!!!!fin de la escritura de las tablas polares!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!!!!!!!!!!!!grafica 5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit= 9, file = "datos.plt", status = "unknown")
  write(9,*)'m = "./data_polares"'!indicaciones del archivo de datos
  write(9,*) 'set title "Energia en el tiempo"'
  write(9,*) 'set autoscale'
  write(9,*) 'unset key'
  write(9,*) 'set xlabel "Tiempo"'
  write(9,*) 'set ylabel "Energia total"'
  write(9,*)'plot m using 1:10 with lines title "Energia"'

  !!!!!!!!!!!!!!!!!!!!!!!grafica 5!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit= 10, file = "datos2.plt", status = "unknown")
  write(10,*)'m = "./data_polares"'!indicaciones del archivo de datos
  write(10,*) 'set title "Energia en el tiempo"'
  write(10,*) 'unset key'
  write(10,*) 'set xlabel "Tiempo"'
  write(10,*) 'set ylabel "Energia total"'
  write(10,*)'plot m using 1:10 with lines title "Energia"'
  write(10,*) 'set xzeroaxis lt -1 lw 1.5' ! Dibuja el eje X (Y=0)
  write(10,*) 'set yzeroaxis lt -1 lw 1.5' ! Dibuja el eje Y (X=0)
  
  call system("gnuplot -p datos.plt")!llamar al sistema para usar gnuplot
 !  call system("gnuplot -p datos2.plt")
  call cpu_time(tfin)
  t_cal = tfin-tini
  print '(A, F10.4, A)', " Tiempo de ejecución del CPU: ", t_cal, " se&
       &gundos."
       !!!!!!!!!!!!!!!!!!!!!!! espacios fase !!!!!!!!!!!!!!!!!!!!!!!!!!!!
  close(2)
  open(unit=11, file="fase.plt", status="replace")
  write(11,*) 'm = "./data_polares"'
  write(11,*) 'nt = 4.0'
  write(11,'(a,f12.6)') 'm1 = ', m1
  write(11,'(a,f12.6)') 'm2 = ', m2
  write(11,*) 'mt = m1 + m2'
  ! columnas: 2=r1 3=r2 4=th1 5=th2 6=vr1 7=vr2 8=w1 9=w2
  write(11,*) 'p1(r1,r2,o1,o2,vr2,w1,w2) = mt*r1**2*w1 &
  &- m2*r1*vr2*sin(o1-o2) + m2*r1*r2*w2*cos(o1-o2)'
  write(11,*) 'p2(r1,r2,o1,o2,vr1,w1,w2) = m2*(r2**2*w2 &
  &+ r2*vr1*sin(o1-o2) + r1*r2*w1*cos(o1-o2))'

  write(11,*) 'set multiplot layout 1,2 title "Espacios fase"'
  write(11,*) 'set key top right'
  write(11,*) 'set grid'
  write(11,*) 'set format x "%.2g"; set format y "%.2g"'

  ! --- Grafica 1: momentos angulares (theta vs p_theta) ---
  write(11,*) 'set title "Momentos angulares"'
  write(11,*) 'set xlabel "theta (rad)"; set ylabel "p_theta"'
  write(11,*) 'stats m using 4 nooutput; xmin=STATS_min; xmax=STATS_max'
  write(11,*) 'stats m using 5 nooutput'
  write(11,*) 'xmin = (STATS_min<xmin) ? STATS_min : xmin'
  write(11,*) 'xmax = (STATS_max>xmax) ? STATS_max : xmax'
  write(11,*) 'stats m using (p1($2,$3,$4,$5,$7,$8,$9)) nooutput'
  write(11,*) 'ymin=STATS_min; ymax=STATS_max'
  write(11,*) 'stats m using (p2($2,$3,$4,$5,$6,$8,$9)) nooutput'
  write(11,*) 'ymin = (STATS_min<ymin) ? STATS_min : ymin'
  write(11,*) 'ymax = (STATS_max>ymax) ? STATS_max : ymax'
  write(11,*) 'dx = (xmax-xmin)/nt; dy = (ymax-ymin)/nt'
  write(11,*) 'set xrange [-1:1.2]'
  write(11,*) 'set yrange [-25:25]'
  write(11,*) 'set xtics dx; set ytics dy'
  write(11,*) 'plot m using 4:(p1($2,$3,$4,$5,$7,$8,$9)) with lines &
  &lc rgb "blue" title "lenteja 1", \'
  write(11,*) '     m using 5:(p2($2,$3,$4,$5,$6,$8,$9)) with lines &
  &lc rgb "red" title "lenteja 2"'

  ! --- Grafica 2: radial (r vs vr) de ambas lentejas ---
  write(11,*) 'set title "Radial"'
  write(11,*) 'set xlabel "r"; set ylabel "vr"'
  write(11,*) 'stats m using 2 nooutput; xmin=STATS_min; xmax=STATS_max'
  write(11,*) 'stats m using 3 nooutput'
  write(11,*) 'xmin = (STATS_min<xmin) ? STATS_min : xmin'
  write(11,*) 'xmax = (STATS_max>xmax) ? STATS_max : xmax'
  write(11,*) 'stats m using 6 nooutput; ymin=STATS_min; ymax=STATS_max'
  write(11,*) 'stats m using 7 nooutput'
  write(11,*) 'ymin = (STATS_min<ymin) ? STATS_min : ymin'
  write(11,*) 'ymax = (STATS_max>ymax) ? STATS_max : ymax'
  write(11,*) 'dx = (xmax-xmin)/nt; dy = (ymax-ymin)/nt'
  write(11,*) 'set xtics dx; set ytics dy'
  write(11,*) 'set xrange [0:2.5]'
  write(11,*) 'set yrange [-24:24]'
  write(11,*) 'plot m using 2:6 with lines lc rgb "blue" title "lenteja 1", \'
  write(11,*) '     m using 3:7 with lines lc rgb "red" title "lenteja 2"'

  write(11,*) 'unset multiplot'
  close(11)

  call system("gnuplot -p fase.plt")
contains

  function dv(t, ve) result(dve)
    use constantes
    use varglob
    implicit none
    real(ki), intent(in) :: t
    real(ki), intent(in) :: ve(8)
    real(ki) :: dve(8)
    
    ! Matrices para resolver el sistema de 4x4
    real(ki) :: M(4,4)
    real(ki) :: Bm(4)
    real(ki) :: a(4) ! Vector para guardar las aceleraciones
    
    ! Variables auxiliares para que tus fórmulas sean más fáciles de leer
    real(ki) :: r1, o1, r2, o2, vr1, w1, vr2, w2, a1,mt
    
    !Desempaquetar variables
    r1 = ve(1); o1 = ve(2); r2 = ve(3); o2 = ve(4)
    vr1 = ve(5); w1 = ve(6); vr2 = ve(7); w2 = ve(8)
    a1 = o1 - o2
    mt = m1+m2

    M = 0.0_ki ! Inicializar todo en cero
    a=0.0_ki
    
    !llenado del sistema
    M(1,1) = mt                             ! Coeficiente de r1_doble_punto
    M(1,2) = 0.0_ki                         ! Coeficiente de theta1_doble_punto
    M(1,3) = m2 * cos(a1)                   ! Coeficiente de r2_doble_punto
    m(1,4) = m2*r2*sin(a1)                  ! Coeficiente de theta2_doble_punto
    
    M(2,1) = 0.0_ki                         
    M(2,2) = mt*(r1**2)              
    M(2,3) = -m2*r1*sin(a1)                  
    M(2,4) = m2*r1*r2*cos(a1)  

    M(3,1) = m2*cos(a1)                   
    M(3,2) = -m2*r1*sin(a1)                         
    M(3,3) = m2                   
    M(3,4) = 0.0_ki

    M(4,1) = r2*sin(a1)                    
    M(4,2) = r1*r2*cos(a1)                        
    M(4,3) = 0.0_ki          
    M(4,4) = (r2**2)  

    
    Bm(1) = (m2*w2)*(r2*w2*cos(a1)-2*vr2*sin(a1))+mt*(r1*(w1**2) + g*cos(o1)) -  kr1*(r1-r01)
    Bm(2) = -mt*(2*r1*vr1*w1+g*r1*sin(o1))-m2*r1*w2*(2*vr2*cos(a1)+r2*w2*sin(a1))
    Bm(3) = m2*(r1*(w1**2)*cos(a1) + r2*(w2**2) + 2*vr1*w1*sin(a1) + g*cos(o2)) - kr2*(r2-r02)
Bm(4) = -r2*(2*vr2*w2 + g*sin(o2) - r1*(w1**2)*sin(a1) + 2*vr1*w1*cos(a1))

    !resolver el sistema
    call resolver_sistema(M, Bm, a)

    ! 5. Ensamblar el vector de derivadas de salida
    dve(1) = vr1    ! derivada de r1 es vr1
    dve(2) = w1     ! derivada de o1 es w1
    dve(3) = vr2
    dve(4) = w2
    
    dve(5) = a(1)   ! derivada de vr1 es a_r1
    dve(6) = a(2)   ! derivada de w1 es a_o1
    dve(7) = a(3)
    dve(8) = a(4)

  end function dv


  !subrutina auspiciada por gemini
subroutine resolver_sistema(M_in, B_in, X)
    use constantes
    implicit none
    real(ki), intent(in)  :: M_in(4,4)
    real(ki), intent(in)  :: B_in(4)
    real(ki), intent(out) :: X(4)
    
    real(ki) :: M(4,4), B(4), temp_row(4), temp_b
    real(ki) :: factor, max_val
    integer  :: i, j, k, max_idx
    
    M = M_in
    B = B_in
    
    do k = 1, 3
       ! ---- PIVOTEO PARCIAL: Busca el valor más grande en la columna ----
       max_val = abs(M(k,k))
       max_idx = k
       do i = k + 1, 4
          if (abs(M(i,k)) > max_val) then
             max_val = abs(M(i,k))
             max_idx = i
          endif
       enddo
       ! Intercambia filas si es necesario
       if (max_idx /= k) then
          temp_row = M(k,:)
          M(k,:) = M(max_idx,:)
          M(max_idx,:) = temp_row
          temp_b = B(k)
          B(k) = B(max_idx)
          B(max_idx) = temp_b
       endif
       ! -----------------------------------------------------------------
       
       if (abs(M(k,k)) < 1.0e-12_ki) then
          X = 0.0_ki
          return ! Evita la división si la matriz de plano es singular
       endif
       
       do i = k + 1, 4
          factor = M(i,k) / M(k,k)
          do j = k, 4
             M(i,j) = M(i,j) - factor * M(k,j)
          end do
          B(i) = B(i) - factor * B(k)
       end do
    end do
    
    X(4) = B(4) / M(4,4)
    do i = 3, 1, -1
       X(i) = B(i)
       do j = i + 1, 4
          X(i) = X(i) - M(i,j) * X(j)
       end do
       X(i) = X(i) / M(i,i)
    end do
  end subroutine resolver_sistema

end program pendulo_doble_elastico_v2
