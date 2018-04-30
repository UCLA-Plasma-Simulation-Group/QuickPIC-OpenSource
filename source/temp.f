      call initialization()                  
      
      do i = start3d, nstep3d

         write (erstr,*) '3D step:', i        
         call perr%werrfl0(erstr)
         
         call MPI_WAIT(id(1),istat,ierr)
         call MPI_WAIT(id(3),istat,ierr)
         call qeb%as(0.0)
         do m = 1, sim%sim%nbeams
            call beam(m)%qdp(qeb)
         end do
         tag(1) = ntag()
         call qeb%ag(tag(1),tag(1),id(1))
         tag(1) = ntag()
         call qeb%pcg(tag(1),tag(1),id(2),id(3))    

         do l =  1, sim%sim%nspecies
            tag_spe(l) = ntag()
            call spe(l)%precv(qe0(l),tag_spe(l))
         end do
         tag(2) = ntag()
         call qi%precv(tag(2))
         tag(3) = ntag()
         call cu%precv(tag(3))
         tag(4) = ntag()
         call epw%precv(tag(4))
         tag(5) = ntag()
         call fxy%precv(tag(5))
         tag(6) = ntag()
         call psit%precv(tag(6))
         tag(7) = ntag()
         call bxyz%precv(tag(7))

         call fxy%cb(bexyz,1,(/1,2/),(/1,2/))         
         call psit%cb(bexyz,1,(/1/),(/3/))
         call bxyz%cb(bbxyz,1,(/1,2,3/),(/1,2,3/))


         if (pp%getstageid() == 0) then
            do m = 1, sim%sim%nspecies
               call qe0(m)%as(0.0)           
               call spe(m)%qdp(qe0(m))
               call qi%add(qi,qe0(m))
               call qe0(m)%cb(qep(m),1,(/1/),(/1/))         
            end do
         end if
                  
         do j = 1, nstep2d
            write (erstr,*) '2D step:', j
            call perr%werrfl0(erstr)
            if (j == nstep2d) then
               call MPI_WAIT(id(2),istat,ierr)
            endif
            call qb%cp(qeb,j+1,(/1/),(/1/))
            call qb%fftrk(1)
            call qb%elf(epwb)
            call qe%mult(qi,-1.0)
            do l = 1, sim%sim%nspecies
               call qe0(l)%as(0.0)
               call spe(l)%qdp(qe0(l))                     
               call qe%add(qe,qe0(l))
            end do
            call qe%fftrk(1)
            call qe%pot(psi)
            call psi%grad(fxy)
            call fxy%fftkr(1)
            call psi%fftkr(1)
            do l = 1, sim%sim%nspecies            
               call spe(l)%extpsi(psi,dex)
            end do
            call psi%mult(psi,dex*dex)
            if (allocated(psi3d)) call psi%cb(psi3d,j+1,(/1/),(/1/))
            call cu%div(div_vpot)
            call div_vpot%pot(psit)
            call psit%mult(psit,-dex)
            call psit%fftkr(1)
            call cu%bf(bxyz)
            do l = 1, iter
               call bxyz%mult(epw,(/1,2/),(/2,1/),(/-dex,dex/))
               call bxyz%mult(bxyz,(/3/),(/3/),(/dex/))
               call bxyz%sub(bxyz,epwb,(/1/),(/1/),(/2/))
               call bxyz%add(bxyz,epwb,(/2/),(/2/),(/1/))
               call bxyz%fftkr(2)
               call cu%as(0.0)
               call dcu%as(0.0)
               call amu%as(0.0)
               do m = 1, sim%sim%nspecies
                  call cu0(m)%as(0.0)
                  call dcu0(m)%as(0.0)
                  call amu0(m)%as(0.0)
                  call spe(m)%amjdp(fxy,bxyz,psit,cu0(m),amu0(m),dcu0(m),dex)
                  call cu0(m)%mult(cu0(m),dex)
                  call amu0(m)%mult(amu0(m),dex)
                  call dcu0(m)%mult(dcu0(m),dex)
                  call cu%add(cu,cu0(m))
                  call amu%add(amu,amu0(m))
                  call dcu%add(dcu,dcu0(m))
               end do
               if (l == iter) then
                  do m = 1, sim%sim%nspecies              
                     call reg%add(qe0(m),cu0(m),(/1/),(/1/),(/3/))
                     call reg%fftrk(1)
                     call reg%smooth(reg)
                     call reg%fftkr(1)
                     call qe0(m)%as(reg)
                     call qe0(m)%cb(qep(m),j+1,(/1/),(/1/))                     
                  end do
                  if (allocated(cu3d)) then
                     call cu%cb(cu3d,j+1,(/1,2,3/),(/1,2,3/))
                  end if               
               endif
               call cu%fftrk(1)
               call dcu%fftrk(1)
               call amu%fftrk(3)
               call cu%bfqp(dcu,amu,epw,dex2,dex)
               call cu%div(div_vpot)
               call div_vpot%pot(psit)
               call psit%mult(psit,-dex)
               call psit%fftkr(1)
               call cu%bf(bxyz)
            enddo
            call bxyz%mult(epw,(/1,2/),(/2,1/),(/-dex,dex/))
            call bxyz%mult(bxyz,(/3/),(/3/),(/dex/))
            call bxyz%sub(bxyz,epwb,(/1/),(/1/),(/2/))
            call bxyz%add(bxyz,epwb,(/2/),(/2/),(/1/))
            call bxyz%fftkr(2)
            call fxy%mult(fxy,-1.0)
            call fxy%add(fxy,bxyz,(/1/),(/1/),(/2/))
            call fxy%sub(fxy,bxyz,(/2/),(/2/),(/1/))
            call dcu%mult(dcu,dxi)
            call cu%sub(cu,dcu,(/1,2/),(/1,2/),(/1,2/))
            do m = 1, sim%sim%nspecies              
               call spe(m)%push(fxy,bxyz,psit,dex)
            end do
            call fxy%mult(fxy,dex)
            call fxy%cb(bexyz,j+1,(/1,2/),(/1,2/))
            call psit%cb(bexyz,j+1,(/1/),(/3/))
            call bxyz%mult(bxyz,(/1,2/),(/1,2/),(/dex,dex/))
            call bxyz%cb(bbxyz,j+1,(/1,2,3/),(/1,2,3/))
         enddo
         
         do m = 1, sim%sim%nspecies                       
            call spe(m)%psend(tag_spe(m),id_spe(m))
         end do
         call MPI_WAIT(id(4),istat,ierr)
         call qi%psend(tag(2),id(4))
         call MPI_WAIT(id(5),istat,ierr)
         call cu%psend(tag(3),id(5))
         call MPI_WAIT(id(6),istat,ierr)
         call epw%psend(tag(4),id(6))
         call MPI_WAIT(id(7),istat,ierr)
         call fxy%psend(tag(5),id(7))
         call MPI_WAIT(id(8),istat,ierr)
         call psit%psend(tag(6),id(8))
         call MPI_WAIT(id(9),istat,ierr)
         call bxyz%psend(tag(7),id(9))

         call diagnostic()

         do m = 1, sim%sim%nbeams
            tag_beam(m) = ntag()
            call MPI_WAIT(id_beam(m),istat,ierr)
            call beam(m)%push(bexyz,bbxyz,dex,dxi,tag_beam(m),tag_beam(m),id_beam(m))
         end do
         
         do m = 1, sim%sim%nspecies                       
            call MPI_WAIT(id_spe(m),istat,ierr)
            call spe(m)%renew(pf(m),qe0(m))
         end do
                          
      enddo
