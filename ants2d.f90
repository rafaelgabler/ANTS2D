program formigas2d

! Declarando as variáveis

integer i,j,k,l				
integer n,m				
integer Nf, Nl				
real xmin,xmax,ymin,ymax		
real dx, dy, dt, Massa, Atrito		
real, allocatable :: S(:,:), F(:,:)     
integer, allocatable :: SAUX(:)		
real, allocatable :: V(:,:), X(:,:) 	
real, allocatable :: RAND(:,:)  	
real, allocatable :: RANDFERO(:,:) 
real, allocatable :: FORCA(:,:) 
real, allocatable :: FICA(:) 
real EVAP(1) 
real, allocatable :: KICKS(:,:) 
integer, allocatable :: RANDFEROI(:,:) 
integer npast				
real p,q				
real tempo				

! i,j,k = Variáveis inteiras utilizadas em loops
! n,m =  Número de divisões em x e y
! Nf = n de formigas
! Nl = n de lattices (células)
! xmin, xmax, etc.. = Vértices do plano de movimento das formigas
! dx, dy, dt, Massa, Atrito = Passo espacial, temporal, massa e atrito das formigas
! S(:,:) = n. de formigas por célula 
! F(:,:) = marcador do feromonio por célula
! SAUX(:) = Vetor auxiliar para identificar qual formiga se encontra em cada celula
! V(:,:) e X(:,:) = Velocidade e posição de cada formiga nas 2 direções do espaço
! RAND(:,:) = Matriz randomica para a simulacao do movimento das formigas
! FORCA(:,:) = Força randomica para movimentar as formigas
! RANDFERO(:,:) = Matriz randomica para povoar o tabuleiro com feromonios (inicial)
! RANDFEROI(:,:) = Matriz randomica para povoar o tabuleiro com feromonios (inicial) - inteira
! FICA(:) = Número randômico entre 0 e 1 que avalia as chances da formiga querer ficar numa nova célula 
! KICKS(:,:) = Matriz randomica para kickar as formigas de células "proibidas" pelas leis
! npast = Numero de passos de tempo
! p = probabilidade do feromonio evaporar 
! q = probabilidade da formiga ficar numa nova célula vazia
! tempo = tempo total de simulação

! Definindo as variáveis

n=20
m=20
Nl=n*m
Nf=25
xmin=0
ymin=0
xmax=10
ymax=10
dt=0.01
Massa=0.1
Atrito=0.1
p=0.6
q=0.8
tempo=20
npast=tempo/dt

dx= (xmax-xmin)/(n-1)
dy= (ymax-ymin)/(m-1)

allocate(S(n,m))
allocate(SAUX(n*m))
allocate(F(n,m))
allocate(V(Nf,2))
allocate(X(Nf,2))
allocate(RAND(Nf,2))
allocate(RANDFERO(n,2))
allocate(RANDFEROI(n,2))
allocate(FORCA(Nf,2))
allocate(FICA(Nf))
allocate(KICKS(Nf,2))

! Distribuição inicial das formigas

call randomica(xmin,xmax,RAND(:,1),Nf,1)
call randomica(ymin,ymax,RAND(:,2),Nf,2)

do i=1,Nf
X(i,1)= RAND(i,1)
X(i,2)= RAND(i,2)
end do

! Montando o campo S (numero de formigas por célula) 

do i=1,Nf

do k=1,n
do j=1,m
if(X(i,1).ge.(k-1)*dx.and.X(i,1).le.k*dx) then
if(X(i,2).ge.(j-1)*dy.and.X(i,2).le.j*dy) then
S(k,j)=S(k,j)+1
SAUX(k*j)=i
end if
end if
end do
end do

end do

! Checando se existe mais de uma formiga por célula

10 do k=1,n
   do j=1,m
   if(S(k,j).gt.1) then
   ! Dou um chega pra lá na direção x da última formiga que entrou na célula de meio dx
   X(SAUX(k*j),1)=X(SAUX(k*j),1)+0.5*dx
   go to 10   
   end if
   end do 
   end do

! Distribuindo inicialmente os feromônios em células aleatórias (rastros de formigas anteriores)

call randomica(1.0,20.0,RANDFERO(:,1),n,3)
call randomica(1.0,20.0,RANDFERO(:,2),n,4)

RANDFEROI=RANDFERO

do k=1,n
do j=1,m
do i=1,NF
if(k.eq.RANDFEROI(i,1).and.j.eq.RANDFEROI(i,2))then
F(k,j)= 1.0
end if
end do
end do
end do


! Abrindo o arquivo de dados para escrever as posições das formigas no tempo

open (2,file='trajetorias.plt')
write(2,*) 'Variables="x","y","u","v"'

open (3,file='feromonios.plt')
write(3,*) 'Variables="x","y"'

! Fazendo o loop principal do código

do k=1,npast


write(*,*) 'Passo',k

! Gerando os números randômicos para computar as forças atuantes em cada formiga

 call randomica(-dx,dx,FORCA(:,1),NF,k) ! Direção x (note o intervalo -dx,dx)
 call randomica(-dy,dy,FORCA(:,2),NF,k+1) ! Direção y (note o intervalo -dy,dy)

! Documentando no arquivo de texto o numero da zona de tempo

write(2,*) 'zone t="',k,'"'
write(3,*) 'zone t="',k,'"'

! Resolvendo velocidade e posição de cada formiga nesse passo de tempo

do i=1,Nf 
 
! Primeiro calculamos a velocidade 

 call resvel(V(i,1),dt,FORCA(i,1),Atrito,Massa) ! Direção x
 call resvel(V(i,2),dt,FORCA(i,2),Atrito,Massa) ! Direção y

! Em seguida calculamos a posição
 
 call respos(X(i,1),dt,V(i,1))	! Direção x
 call respos(X(i,2),dt,V(i,2))	! Direção y

! Evitando que as formigas escapem do tabuleiro

if(X(i,1).le.xmin) then
X(i,1)=X(i,1)+0.1*dx
end if
if(X(i,1).ge.xmax) then
X(i,1)=X(i,1)-0.1*dx
end if
if(X(i,2).le.ymin) then
X(i,2)=X(i,2)+0.1*dy
end if
if(X(i,2).ge.ymax) then
X(i,2)=X(i,2)-0.1*dy
end if

end do

! Aqui devemos começar o processo de linkar os rastros de feromonios com o movimento delas

! Atualizando o campo S (numero de formigas por célula) 

20 S=0.0

  do i=1,NF
    do l=1,n
     do j=1,m
       if(X(i,1).ge.(l-1)*dx.and.X(i,1).le.l*dx) then
        if(X(i,2).ge.(j-1)*dy.and.X(i,2).le.j*dy) then
         S(l,j)=S(l,j)+1
         SAUX(l*j)=i
        end if
       end if


! Checando se existe mais de uma formiga por célula

         if(S(l,j).gt.1) then

   ! Se existir, dou um chega pra lá na direção x da última formiga que entrou na célula 
  
         X(SAUX(l*j),1)=X(SAUX(l*j),1)+ FORCA(i,1)*0.1
         X(SAUX(l*j),2)=X(SAUX(l*j),2)+ FORCA(i,2)*0.1
               go to 20   
         end if 
       end do
     end do
    end do

!  Vamos agora verificar se a formiga vai querer ficar na nova célula
!  Isso vai depender da existência de feromônios dentro dessa célula!
   
! Mas antes de tudo vamos pegar um número aleatório entre 0 e 1 para cada formiga

call randomica(0.0,1.0,FICA,Nf,npast+2)
  
! Agora vamos checar para cada formiga se essa nova célula ocupada possui feromonio 
    do l=1,n ! Fazendo uma varredura em todas as células do tabuleiro
     do j=1,m ! Fazendo uma varredura em todas as células do tabuleiro
       if(S(l,j).eq.1) then ! Se aquela célula estiver ocupada por uma formiga
        if(F(l,j).eq.0) then ! E se NÃO existir feromonio nela     
! Então ela tem uma chance "s" de ficar, sendo s < q, onde q é a chance dela ficar se tiver feromonio na célula. Para fazer isso, vamos considerar que cada formiga terá um número randômico entre 0 e 1 vinculado a ela. Esses números estão armazenados no vetor FICA(i), em que i = 1,Nf. Se esse número for menor que q (chance "q" de ocorrer) a formiga fica na célula, se for maior que q (chance "s=1-q" de ocorrer) a formiga vai ter que vazar da célula (aqui darei simplesmente um chega pra lá nela)          
	  if(FICA(SAUX(l*j)).gt.q) then
	    	X(SAUX(l*j),1)=X(SAUX(l*j),1)+ FORCA(SAUX(l*j),1)*0.1
		X(SAUX(l*j),2)=X(SAUX(l*j),2)+ FORCA(SAUX(l*j),2)*0.1

          else    

		X(SAUX(l*j),1)=X(SAUX(l*j),1)
		X(SAUX(l*j),2)=X(SAUX(l*j),2)
	  end if

          end if
          end if
       end do
     end do

! Agora precisamos deixar um rastro de feromônio nas células ocupadas por formigas

do l=1,n
do j=1,m
   if(S(l,j).eq.1) then ! Se tiver uma formiga na célula
   F(l,j)=1.0 ! Então colocamos um feromonio ali
   else ! Caso contrário vamos ver o feromonio vai ou nao evaporar
   if(F(l,j).eq.1.0) then ! Aqui é o caso de termos uma célula vazia com feromônio   
   call randomica(0.0,1.0,EVAP,1,npast+2*k)  ! Pedimos um numero randomico entre 0 e 1 
   if(EVAP(1).lt.p) then ! Se esse numero for menor que p (chance "p" de ocorrer)
   F(l,j)=0.0 ! O feromonio evapora
   end if
   end if
   end if
end do
end do
  
! Escrevendo a posição e a velocidade das formigas num arquivo de dados

do i=1,Nf
write(2,'(F12.4,F12.4,F12.4,F12.4)') X(i,1),X(i,2),V(i,1),V(i,2)
end do

! Escrevendo a posição dos feromônios num arquivo de dados
do l=1,n ! Fazemos uma varredura em x
do j=1,m ! E em y
if(F(l,j).eq.1.0) then ! E se tiver um feromonio naquela célula
write(3,'(F12.4,F12.4)') l*dx+0.5,j*dy+0.5  ! Imprimimos as coordenadas do centro da célula
end if
end do
end do

! Finalizando o "do" vinculado ao passo de tempo (k)
end do

end

!*****************************************************************************************!
!                        Subrotinas utilizadas no código                                  !
!*****************************************************************************************!

! Resolve a velocidade de cada formiga

subroutine resvel(a,b,c,d,e)
real a                      ! velocidade
real b			    ! passo de tempo
real c                      ! forca randomica
real d			    ! atrito
real e			    ! massa
real k1,k2,k3,k4            ! variaveis internas utilizadas para o runge-kutta de 4 ordem

k1 = b*(c-d*a)/e
k2 = b*(c-d*(a+(k1*0.5)))/e
k3 = b*(c-d*(a+(k2*0.5)))/e
k4 = b*(c-d*(a+k3))/e

a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
end subroutine resvel

! Resolve a posição de cada formiga

subroutine respos(a,b,c)
real a                      ! posicao
real b			    ! passo de tempo
real c                      ! componente de velocidade em questao
real k1,k2,k3,k4            ! variaveis internas utilizadas para o runge-kutta de 4 ordem

k1=b*(c)
k2=b*((0.5*k1)+c)
k3=b*((0.5*k2)+c)
k4=b*((k3)+c)
a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
end subroutine respos

! Subrotina de geração dos números randômicos

subroutine randomica(a,b,c,n,d)
real a,b                    ! a,b = range do numero randomico
integer n, m
real c(n)                   ! c = sequencia randomica gerada
integer d,i,e
integer f(8)
integer, allocatable :: seed(:)

 call random_seed(size = m)
allocate (seed(m))

 CALL DATE_AND_TIME(values=f)
 CALL SYSTEM_CLOCK(count=e)

do i = 1,m
seed(i) =  47*d + f(8)*i*d*12000 + e*(3*d+i)
end do

 call random_seed(put = seed)

 call random_number(c)

 c = a+(b-a)*c

end subroutine randomica
