*************************************************************************
*     This program is to identify the best alignment of two protein 
*     structures to give the best TM-score. By default, TM-score is 
*     normalized by the second protein. The program can be freely 
*     copied or modified or redistributed.
*     (For comments, please email to: yzhang@ku.edu)
*     
*     Reference:
*     Yang Zhang, Jeffrey Skolnick, Nucl. Acid Res. 2005 33: 2303-9
*
************************ updating history *******************************
*     2005/06/01: A small bug of two-point superposition was fixed.
*     2005/10/19: the program was reformed so that the alignment
*                 results are not dependent on the specific compilers.
*     2006/06/20: select 'A' if there is altLoc when reading PDB file.
*     2007/02/27: rotation matrix from Chain-1 to Chain-2 is added.
*     2007/04/18: add options with TM-score normalized by average
*                 length, shorter length, or longer length of two 
*                 structures.
*     2007/05/23: add additional output file 'TM.sup_all' for showing
*                 all atoms while 'TM.sup' is only for aligned atoms
*     2007/09/19: add a new feature alignment to deal with the problem
*                 of aligning fractional structures (e.g. protein
*                 interfaces).
*     2007/10/16: A bug for irregular bond-length models was fixed
*     
*     2008/07/21: Modified options to allow file-based output of various
*                 results: aligned residues, scores, xform matrix, etc.
*
*     Version patched by SBGrid at Harvard Medical School
*     Web:    http://www.sbgrid.org
*     Source: $HeadURL: https://sbgrid.org:4443/svn/ijstokes/trunk/projects/tmalign/src/TMalign.f $
*     Id:     $Id: TMalign.f 280 2009-08-13 15:40:01Z ijstokes $ 
*************************************************************************
      
      program compares
      PARAMETER(nmax=5000)
      PARAMETER(nmax2=10000)

      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/length/nseq1,nseq2
      common/d0/d0,anseq
      common/d0min/d0_min
      common/d00/d00,d002
      
      character atoms(10000)*100
      integer   res_end(nmax)
      integer   atom_idx, res_idx, cur_res_idx, t_idx
      integer   pdb_len, verbose
      integer   save_all_res,save_xform,save_score,save_align
      integer   print_all_res,print_xform,print_score,print_align,print_dis

      character*120 fnam,pdb(100),outname,pdbouttmp
      character*120 pdbname,scorename,xformname,alignname
      character*20  label
      character*13  pdbfragout
      character*3 aa(-1:20),aanam,ss1(nmax),ss2(nmax)
      character*120 s,du
      character*200 outnameall_tmp,outnameall
      character seq1(0:nmax),seq2(0:nmax)
      character aseq1(nmax2),aseq2(nmax2),aseq3(nmax2)

      dimension xx(nmax),yy(nmax),zz(nmax)
      dimension m1(nmax),m2(nmax)
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/init/invmap_i(nmax)

      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d8/d8
      common/initial4/mm1(nmax),mm2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      double precision origx,origy,origz
      double precision newx,newy,newz
      character*30 pre
      character*26 post
      data w /nmax*1.0/
ccc   

      data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     &     'PRO','MET','ASP','ASN','LEU',
     &     'LYS','GLU','GLN','ARG',
     &     'HIS','PHE','TYR','TRP','CYX'/
      character*1 slc(-1:20)
      data slc/'X','G','A','S','C','V','T','I',
     &     'P','M','D','N','L','K','E','Q','R',
     &     'H','F','Y','W','C'/

      call getarg(1,fnam)
      if(fnam.eq.' '.or.fnam.eq.'?'.or.fnam.eq.'-h')then
         write(*,*)'tmalign struct1.pdb struct2.pdb [options]'
         write(*,*)
         write(*,*)'-h           help'
         write(*,*)'-q           quiet (supress screen output)'
         write(*,*)
         write(*,*)'-r <filname> save aligned residues to PDB file'
         write(*,*)'-x <filname> save transform matrix'
         write(*,*)'-s <filname> save scores'
         write(*,*)'-a <filname> save alignment strings'

         write(*,*)'-R output aligned residues to PDB file'
         write(*,*)'-X output transform matrix'
         write(*,*)'-S output scores'
         write(*,*)'-A output alignment strings'
         write(*,*)'-D output residue pair distances'
         
         write(*,*)
         write(*,*)'-l <label>   label appended to score output'
         write(*,*)
         write(*,*)'-ave         normalize score by avg chain length'
         write(*,*)'-short       normalize score by short chain length'
         write(*,*)'-long        normalize score by long chain length'
         write(*,*)'-L <length>  normalize score by arbitrary length'
         write(*,*)
         write(*,*)'Brief instruction for running TM-align program:'
         write(*,*)'(For detail: Zhang & Skolnick, Nucl. Acid Res.',
     &     '2005 33, 2303)'
         write(*,*)
         goto 9999
      endif

******* options ----------->
      save_all_res = 0          !whether or not to output all residues with alignment
      save_xform = 0
      save_score = 0
      save_align = 0

      print_all_res = 0          !whether or not to output all residues with alignment
      print_xform = 0
      print_score = 0
      print_align = 0
      print_dis   = 0
      
      verbose = 1
      label = ''
      m_out=-1                  !decided output
      m_fix=-1                  !fixed length-scale only for output
      m_ave=-1                  !using average length
      m_d0_min=-1               !diminum d0 for search
      m_d0=-1                   !given d0 for both search and output
      narg=iargc()
      i=0
      j=0
      atom_idx=0
      res_idx=-1
      cur_res_idx=0
 115  continue
      i=i+1
      call getarg(i,fnam)
      if (fnam.eq.'-r') then
         save_all_res=1
         i=i+1
         call getarg(i,pdbname)
      elseif (fnam.eq.'-x') then
         save_xform=1
         i=i+1
         call getarg(i,xformname)
      elseif (fnam.eq.'-s') then
         save_score=1
         i=i+1
         call getarg(i,scorename)
      elseif (fnam.eq.'-a') then
         save_align=1
         i=i+1
         call getarg(i,alignname)
      elseif (fnam.eq.'-R') then
         print_all_res=1
      elseif (fnam.eq.'-X') then
         print_xform=1
      elseif (fnam.eq.'-S') then
         print_score=1
      elseif (fnam.eq.'-A') then
         print_align=1
      elseif (fnam.eq.'-D') then
         print_dis=1
      elseif (fnam.eq.'-l') then
         i=i+1
         call getarg(i,label)
      elseif (fnam.eq.'-L') then  !change both L_all and d0
         m_fix=1
         i=i+1
         call getarg(i,fnam)
         read(fnam,*)L_fix
      elseif (fnam.eq.'-ave') then ! this will change the superposed output but not the alignment
         m_ave=1
      elseif (fnam.eq.'-short') then
         m_ave=2
      elseif (fnam.eq.'-long') then
         m_ave=3
      elseif (fnam.eq.'-q') then
         verbose = 0
      else
         j=j+1
         pdb(j)=fnam
      endif
      
      if(i.lt.narg)goto 115
      
ccccccccc read data from first CA file:
      open(unit=10,file=pdb(1),status='old')
      i=0
      do while (.true.)
         read(10,9001,end=1010) s

         if(i.gt.0.and.s(1:3).eq.'TER')goto 1010

         if(s(1:3).eq.'ATO')then
            atom_idx = atom_idx + 1
            atoms(atom_idx) = s

            read(unit=s(23:26),fmt='(I4)') cur_res_idx
            if(cur_res_idx .ne. res_idx .or. s(27:27) .ne. ' ') then
                res_end(i) = atom_idx - 1
                res_idx = cur_res_idx
            endif

            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or
     &           .s(13:16).eq.'  CA')then
               if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
                  i=i+1
                  read(s,9000)du,aanam,du,mm1(i),du,
     $                 xa(1,i,0),xa(2,i,0),xa(3,i,0)
                  do j=-1,20
                     if(aanam.eq.aa(j))seq1(i)=slc(j)
                  enddo
                  ss1(i)=aanam
                  if(i.ge.nmax)goto 1010
               endif
            endif
         endif
      enddo
 1010 continue
      res_end(i) = atom_idx
 9000 format(A17,A3,A2,i4,A4,3F8.3)
 9001 format(A100)
      close(10)
      nseq1=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ccccccccc read data from the second CA file:
      open(unit=11,file=pdb(2),status='old')
      i=0
      do while (.true.)
         read(11,9001,end=1011) s
         if(i.gt.0.and.s(1:3).eq.'TER')goto 1011
         if(s(1:3).eq.'ATO')then

            if(s(13:16).eq.'CA  '.or.s(13:16).eq.' CA '.or.
     &           s(13:16).eq.'  CA')then
c              we are now at the terminal C-alpha
               if(s(17:17).eq.' '.or.s(17:17).eq.'A')then
                  i=i+1
                  read(s,9000)du,aanam,du,mm2(i),du,
     $                 xa(1,i,1),xa(2,i,1),xa(3,i,1)
                  do j=-1,20
                     if(aanam.eq.aa(j))seq2(i)=slc(j)
                  enddo
                  ss2(i)=aanam
                  if(i.ge.nmax)goto 1011
               endif
            endif
         endif
      enddo
 1011 continue
      close(11)
      nseq2=i
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

*!!!  Scale of TM-score in search is based on the smaller protein --------->
      d0_min=0.5
      if(m_d0_min.eq.1)then
         d0_min=d0_min_input    !for search
      endif
      anseq_min=min(nseq1,nseq2)
      anseq=anseq_min           !length for defining TMscore in search
      d8=1.5*anseq_min**0.3+3.5 !remove pairs with dis>d8 during search & final
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix
      d00=d0                    !for quickly calculate TM-score in searching
      if(d00.gt.8)d00=8
      if(d00.lt.4.5)d00=4.5
      d002=d00**2
      nseq=max(nseq1,nseq2)
      do i=1,nseq
         n1(i)=i
         n2(i)=i
      enddo
      
***** do alignment **************************
      CALL super_align          !to find invmap(j)
      
************************************************************
***   resuperpose to find residues of dis<d8 ------------------------>
      n_al=0
      do j=1,nseq2
         if(invmap0(j).gt.0)then
            i=invmap0(j)
            n_al=n_al+1
            xtm1(n_al)=xa(1,i,0)
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
            m1(n_al)=i          !for recording residue order
            m2(n_al)=j
         endif
      enddo
      d0_input=d0
      call TMscore8(d0_input,n_al,xtm1,ytm1,ztm1,n1,n_al,
     &     xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !TM-score with dis<d8 only

*!!!  Output TM-score is based on the second protein------------------>
      d0_min=0.5                !for output
      anseq=nseq2               !length for defining final TMscore
      if(m_ave.eq.1)anseq=(nseq1+nseq2)/2.0 !<L>
      if(m_ave.eq.2)anseq=min(nseq1,nseq2)
      if(m_ave.eq.3)anseq=max(nseq1,nseq2)
      if(anseq.lt.anseq_min)anseq=anseq_min
      if(m_fix.eq.1)anseq=L_fix !input length
      if(anseq.gt.15)then
         d0=1.24*(anseq-15)**(1.0/3.0)-1.8 !scale for defining TM-score
      else
         d0=d0_min
      endif
      if(d0.lt.d0_min)d0=d0_min
      if(m_d0.eq.1)d0=d0_fix
      
***   remove dis>d8 in normal TM-score calculation for final report----->
      j=0
      n_eq=0
      do i=1,n_al
         dis2=sqrt((xtm1(i)-xtm2(i))**2+(ytm1(i)-ytm2(i))**2+
     &        (ztm1(i)-ztm2(i))**2)
         if(dis2.le.d8)then
            j=j+1
            xtm1(j)=xtm1(i)
            ytm1(j)=ytm1(i)
            ztm1(j)=ztm1(i)
            xtm2(j)=xtm2(i)
            ytm2(j)=ytm2(i)
            ztm2(j)=ztm2(i)
            m1(j)=m1(i)
            m2(j)=m2(i)
            if(ss1(m1(i)).eq.ss2(m2(i)))then
               n_eq=n_eq+1
            endif
         endif
      enddo
      seq_id=float(n_eq)/(n_al+0.00000001)
      n8_al=j
      d0_input=d0
      call TMscore(d0_input,n8_al,xtm1,ytm1,ztm1,n1,n8_al,
     &     xtm2,ytm2,ztm2,n2,TM8,Rcomm,Lcomm) !normal TMscore
      rmsd8_al=Rcomm
      TM8=TM8*n8_al/anseq       !TM-score after cutoff
      
********* for output summary ******************************
      if (verbose .eq. 1) then
          write(*,101)nseq1,pdb(1)
          write(*,102)nseq2,pdb(2)
          write(*,108),int(anseq)
          write(*,*)
          write(*,103)n8_al,rmsd8_al,TM8,seq_id      
          write(*,*)
      endif

      if (save_score) then
          open(unit=15,file=scorename,status='replace')
          write(unit=15,fmt=107)n8_al,rmsd8_al,TM8,seq_id,nseq1,nseq2,
     &    label
          close(15)
      endif
      
      if (print_score) then
          write(*,fmt=107)n8_al,rmsd8_al,TM8,seq_id,nseq1,nseq2,label
      endif
     
 101  format('Chain 1 (Size=',I4,'):',A100)
 102  format('Chain 2 (Size=',I4,'):',A100)
 108  format('TM-score is normalized by ',I4)
 103  format('Aligned length=',I4,', RMSD=',f6.2,
     &     ', TM-score=',f7.5,', ID=',f5.3)
 107  format(I4,' ',f6.2,' ',f7.5,' ',f5.3,' ',I4,' ',I4,' ',A20)
 
********* extract rotation matrix ------------>
      L=0
      do i=1,n8_al
         k=m1(i)
         L=L+1
         r_1(1,L)=xa(1,k,0)
         r_1(2,L)=xa(2,k,0)
         r_1(3,L)=xa(3,k,0)
         r_2(1,L)=xtm1(i)
         r_2(2,L)=ytm1(i)
         r_2(3,L)=ztm1(i)
      enddo
      
      
      if(L.gt.3)then
         call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
         armsd=dsqrt(rms/L)
         if (verbose .eq. 1) then
           write(*,*)'-------- rotation matrix to rotate Chain-1 to ',
     &        'Chain-2 ------'
           write(*,*)'i          t(i)         u(i,1)         u(i,2) ',
     &        '        u(i,3)'
           do i=1,3
              write(*,204)i,t(i),u(i,1),u(i,2),u(i,3)
           enddo
           write(*,*)
         endif

c PDBROT
         if (save_xform) then
            open(unit=13,file=xformname,status='replace')
            write(13,205)t(1),' ',t(2),' ',t(3)
            do i=1,3
                write(13,205)u(i,1),' ',u(i,2),' ',u(i,3)
            enddo
            close(13)
         endif
         if (print_xform) then
            write(*,fmt=205)t(1),' ',t(2),' ',t(3)
            do i=1,3
                write(*,fmt=205)u(i,1),' ',u(i,2),' ',u(i,3)
            enddo
         endif
      endif
 204  format(I2,f18.10,f15.10,f15.10,f15.10)
 205  format(f15.10,A1,f15.10,A1,f15.10)

************  output aligned sequences **************************
      if (save_all_res) then
          open(unit=12,file=pdbname,status='replace')
      endif
      
      ii=0
      i1_old=1
      i2_old=1
      
      do i=1,n8_al
         do j=i1_old,m1(i)-1
            ii=ii+1
            aseq1(ii)=seq1(j)
            aseq2(ii)='-'
            aseq3(ii)=' '
         enddo
         do j=i2_old,m2(i)-1
            ii=ii+1
            aseq1(ii)='-'
            aseq2(ii)=seq2(j)
            aseq3(ii)=' '
         enddo
         ii=ii+1
         aseq1(ii)=seq1(m1(i))
         aseq2(ii)=seq2(m2(i))
         dis2=sqrt((xtm1(i)-xtm2(i))**2+
     &     (ytm1(i)-ytm2(i))**2+(ztm1(i)-ztm2(i))**2)
         if(dis2.le.5)then
           aseq3(ii)=':'
c           write(*,*) i, ii, m1(i), m2(i), res_end(m1(i)-1)+1,
c     &                res_end(m1(i)), aseq1(ii), aseq2(ii)
         else
           aseq3(ii)='.'
         endif
         if (print_dis.eq.1) then
           write(*,*) dis2
         endif
         if (save_all_res.or.print_all_res) then
           do t_idx = res_end(m1(i)-1)+1, res_end(m1(i))
              read(atoms(t_idx),9010) pre, origx, origy, origz, post
              newx=t(1)+u(1,1)*origx+u(1,2)*origy+u(1,3)*origz
              newy=t(2)+u(2,1)*origx+u(2,2)*origy+u(2,3)*origz
              newz=t(3)+u(3,1)*origx+u(3,2)*origy+u(3,3)*origz

              if (save_all_res) then
                 write(unit=12,fmt=9010) pre, newx, newy, newz, post
              endif
              
              if (print_all_res) then
                 write(*,fmt=9010) pre, newx, newy, newz, post
              endif
              
           enddo
         endif
         i1_old=m1(i)+1
         i2_old=m2(i)+1
      enddo
      
      if (save_all_res) then
         close(12)
      endif
         
      do i=i1_old,nseq1
         ii=ii+1
         aseq1(ii)=seq1(i)
         aseq2(ii)='-'
         aseq3(ii)=' '
      enddo
      
      do i=i2_old,nseq2
         ii=ii+1
         aseq1(ii)='-'
         aseq2(ii)=seq2(i)
         aseq3(ii)=' '
      enddo
      
      if (verbose .eq. 1) then
          write(*,50)
          write(*,10)(aseq1(i),i=1,ii)
          write(*,10)(aseq3(i),i=1,ii)
          write(*,10)(aseq2(i),i=1,ii)
          write(*,*)
      endif
      
      if (save_align) then
          open(unit=14,file=alignname,status='replace')
          write(unit=14,fmt=10)(aseq1(i),i=1,ii)
          write(unit=14,fmt=10)(aseq3(i),i=1,ii)
          write(unit=14,fmt=10)(aseq2(i),i=1,ii)
          close(14)
      endif
      
      if (print_align) then
          write(*,fmt=10)(aseq1(i),i=1,ii)
          write(*,fmt=10)(aseq3(i),i=1,ii)
          write(*,fmt=10)(aseq2(i),i=1,ii)
      endif
      
 10   format(10000A1)
 50   format('(":" marks residue pairs with a distance < 5.0 Angstrom)')
 9010 format(A30,3F8.3,A26)

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 9999 END

***********************************************************************
***********************************************************************
*     Structure superposition
***********************************************************************
***********************************************************************
***********************************************************************
      SUBROUTINE super_align
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      dimension gapp(100)

      TMmax=0
      n_gapp=2
      gapp(1)=-0.6
      gapp(2)=0

c      n_gapp=11
c      do i=1,n_gapp
c         gapp(i)=-(n_gapp-i)
c      enddo

*11111111111111111111111111111111111111111111111111111111
*     get initial alignment from gapless threading
**********************************************************
      call get_initial          !gapless threading
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, matrix score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

*****************************************************************
*       initerative alignment, for different gap_open:
*****************************************************************
      DO 111 i_gapp=1,n_gapp !different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 222 id=1,30         !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 33
            endif
            TM_old=TM
 222     continue
 33      continue
 111  continue

*222222222222222222222222222222222222222222222222222222222
*     get initial alignment from secondary structure alignment
**********************************************************
      call get_initial2         !DP for secondary structure
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 1111 i_gapp=1,n_gapp !different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 2222 id=1,30 !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
c     write(*,21)gap_open,rmsd_al,n_al,TM
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 333
            endif
            TM_old=TM
 2222    continue
 333     continue
 1111 continue
      
*333333333333333333333333333333333333333333333333333333333333
*     get initial alignment from invmap0+SS
*************************************************************
      call get_initial3         !invmap0+SS
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 1110 i_gapp=1,n_gapp !different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 2220 id=1,30 !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)

            call get_score      !calculate TM-score, score(i,j)
c     write(*,21)gap_open,rmsd_al,n_al,TM
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
            if(id.gt.1)then
               diff=abs(TM-TM_old)
               if(diff.lt.0.000001)goto 330
            endif
            TM_old=TM
 2220    continue
 330     continue
 1110 continue

*444444444444444444444444444444444444444444444444444444444
*     get initial alignment of pieces from gapless threading
**********************************************************
      call get_initial4         !gapless threading
      do i=1,nseq2
         invmap(i)=invmap_i(i)  !with highest zcore
      enddo
      call get_score            !TM, matrix score(i,j)
      if(TM.gt.TMmax)then
         TMmax=TM
         do j=1,nseq2
            invmap0(j)=invmap(j)
         enddo
      endif

*****************************************************************
*     initerative alignment, for different gap_open:
*****************************************************************
      DO 44 i_gapp=2,n_gapp    !different gap panalties
         GAP_OPEN=gapp(i_gapp)  !gap panalty
         do 444 id=1,2          !maximum interation is 200
            call DP(NSEQ1,NSEQ2) !produce alignment invmap(j)
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
            
            call get_score      !calculate TM-score, score(i,j)
c     record the best alignment in whole search ---------->
            if(TM.gt.TMmax)then
               TMmax=TM
               do j=1,nseq2
                  invmap0(j)=invmap(j)
               enddo
            endif
 444     continue
 44   continue

c^^^^^^^^^^^^^^^ best alignment invmap0(j) found ^^^^^^^^^^^^^^^^^^
      RETURN
      END

**************************************************************
*     get initial alignment invmap0(i) from gapless threading
**************************************************************
      subroutine get_initial
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)

      aL=min(nseq1,nseq2)
      idel=aL/2.5               !minimum size of considered fragment
      if(idel.le.5)idel=5
      n1=-nseq2+idel
      n2=nseq1-idel
      GL_max=0
      do ishift=n1,n2
         L=0
         do j=1,nseq2
            i=j+ishift
            if(i.ge.1.and.i.le.nseq1)then
               L=L+1
               invmap(j)=i
            else
               invmap(j)=-1
            endif
         enddo
         if(L.ge.idel)then
            call get_GL(GL)
            if(GL.gt.GL_max)then
               GL_max=GL
               do i=1,nseq2
                  invmap_i(i)=invmap(i)
               enddo
            endif
         endif
      enddo

      return
      end

**************************************************************
*     get initial alignment invmap0(i) from secondary structure
**************************************************************
      subroutine get_initial2
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/sec/isec(nmax),jsec(nmax)
      common/init/invmap_i(nmax)

********** assign secondary structures ***************
c     1->coil, 2->helix, 3->turn, 4->strand
      do i=1,nseq1
         isec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq1)then
            dis13=diszy(0,j1,j3)
            dis14=diszy(0,j1,j4)
            dis15=diszy(0,j1,j5)
            dis24=diszy(0,j2,j4)
            dis25=diszy(0,j2,j5)
            dis35=diszy(0,j3,j5)
            isec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo
      do i=1,nseq2
         jsec(i)=1
         j1=i-2
         j2=i-1
         j3=i
         j4=i+1
         j5=i+2
         if(j1.ge.1.and.j5.le.nseq2)then
            dis13=diszy(1,j1,j3)
            dis14=diszy(1,j1,j4)
            dis15=diszy(1,j1,j5)
            dis24=diszy(1,j2,j4)
            dis25=diszy(1,j2,j5)
            dis35=diszy(1,j3,j5)
            jsec(i)=make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
         endif
      enddo
      call smooth               !smooth the assignment

********** score matrix **************************
      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=1
            else
               score(i,j)=0
            endif
         enddo
      enddo

********** find initial alignment: invmap(j) ************
      gap_open=-1.0             !should be -1
      call DP(NSEQ1,NSEQ2)      !produce alignment invmap(j)
      do i=1,nseq2
         invmap_i(i)=invmap(i)
      enddo

*^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

**************************************************************
*     get initial alignment invmap0(i) from secondary structure 
*     and previous alignments
**************************************************************
      subroutine get_initial3
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/sec/isec(nmax),jsec(nmax)
      common/init/invmap_i(nmax)

********** score matrix **************************
      do i=1,nseq2
         invmap(i)=invmap0(i)
      enddo
      call get_score1           !get score(i,j) using RMSD martix
      do i=1,nseq1
         do j=1,nseq2
            if(isec(i).eq.jsec(j))then
               score(i,j)=0.5+score(i,j)
            else
               score(i,j)=score(i,j)
            endif
         enddo
      enddo

********** find initial alignment: invmap(j) ************
      gap_open=-1.0             !should be -1
      call DP(NSEQ1,NSEQ2)      !produce alignment invmap(j)
      do i=1,nseq2
         invmap_i(i)=invmap(i)
      enddo

*^^^^^^^^^^^^ initial alignment done ^^^^^^^^^^^^^^^^^^^^^^
      return
      end

**************************************************************
*     get initial alignment invmap0(i) from fragment gapless threading
**************************************************************
      subroutine get_initial4
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/length/nseq1,nseq2
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/alignrst/invmap0(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/TM/TM,TMmax
      common/init/invmap_i(nmax)
      common/initial4/mm1(nmax),mm2(nmax)
      logical contin

      dimension ifr2(2,nmax,nmax),Lfr2(2,nmax),Lfr_max2(2),i_fr2(2)
      dimension ifr(nmax)
      dimension mm(2,nmax)
      
      fra_min=4                 !>=4,minimum fragment for search
      fra_min1=fra_min-1        !cutoff for shift, save time
      dcu0=3.85
c      dcu_min=3.65

ccc   Find the smallest continuous fragments -------->
      do i=1,nseq1
         mm(1,i)=mm1(i)
      enddo
      do i=1,nseq2
         mm(2,i)=mm2(i)
      enddo
      do k=1,2
         dcu=dcu0
         if(k.eq.1)then
            nseq0=nseq1
            r_min=nseq1/3.0     !minimum fragment, in case too small protein
         else
            nseq0=nseq2
            r_min=nseq2/3.0     !minimum fragment, in case too small protein
         endif
         if(r_min.gt.fra_min)r_min=fra_min
 20      nfr=1                  !number of fragments
         j=1                    !number of residue at nf-fragment
         ifr2(k,nfr,j)=1        !what residue
         Lfr2(k,nfr)=j          !length of the fragment
         do i=2,nseq0
            dis=diszy(k-1,i-1,i)
            contin=.false.
            if(dcu.gt.dcu0)then
               if(dis.lt.dcu)then
c                  if(dis.gt.dcu_min)then
                     contin=.true.
c                  endif
               endif
            elseif(mm(k,i).eq.(mm(k,i-1)+1))then
               if(dis.lt.dcu)then
c                  if(dis.gt.dcu_min)then
                     contin=.true.
c                  endif
               endif
            endif
            if(contin)then
               j=j+1
               ifr2(k,nfr,j)=i
               Lfr2(k,nfr)=j
            else
               nfr=nfr+1
               j=1
               ifr2(k,nfr,j)=i
               Lfr2(k,nfr)=j
            endif
         enddo
         Lfr_max=0
         i_fr2(k)=1             !ID of the maximum piece
         do i=1,nfr
            if(Lfr_max.lt.Lfr2(k,i))then
               Lfr_max=Lfr2(k,i)
               i_fr2(k)=i
            endif
         enddo
         if(Lfr_max.lt.r_min)then
            dcu=1.1*dcu
            goto 20
         endif
      enddo
c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      
ccc   select what piece will be used (this may araise ansysmetry, but
ccc   only when L1=L2 and Lfr1=Lfr2 and L1 ne Lfr1
ccc   if L1=Lfr1 and L2=Lfr2 (normal proteins), it will be the same as initial1
      mark=1
      if(Lfr2(1,i_fr2(1)).lt.Lfr2(2,i_fr2(2)))then
         mark=1
      elseif(Lfr2(1,i_fr2(1)).gt.Lfr2(2,i_fr2(2)))then
         mark=2
      else                      !Lfr1=Lfr2
         if(nseq1.le.nseq2)then
            mark=1
         else
            mark=2
         endif
      endif
ccc   
      L_fr=Lfr2(mark,i_fr2(mark))
      do i=1,L_fr
         ifr(i)=ifr2(mark,i_fr2(mark),i)
      enddo
ccc   
      if(mark.eq.1)then         !non-redundant to get_initial1
         nseq0=nseq1
      else
         nseq0=nseq2
      endif
      if(L_fr.eq.nseq0)then
         n1=int(nseq0*0.1)      !0
         n2=int(nseq0*0.89)     !2
         j=0
         do i=n1,n2
            j=j+1
            ifr(j)=ifr(n1+j)
         enddo
         L_fr=j
      endif
      
ccc   get initial ------------->
      if(mark.eq.1)then    !nseq1 as the smallest one
         nseq1_=L_fr
         aL=min(nseq1_,nseq2)
         idel=aL/2.5            !minimum size of considered fragment
         if(idel.le.fra_min1)idel=fra_min1
         n1=-nseq2+idel         !shift1
         n2=nseq1_-idel         !shift2
         GL_max=0
         do ishift=n1,n2
            L=0
            do j=1,nseq2
               i=j+ishift
               if(i.ge.1.and.i.le.nseq1_)then
                  L=L+1
                  invmap(j)=ifr(i)
               else
                  invmap(j)=-1
               endif
            enddo
            if(L.ge.idel)then
               call get_GL(GL)
               if(GL.gt.GL_max)then
                  GL_max=GL
                  do i=1,nseq2
                     invmap_i(i)=invmap(i)
                  enddo
               endif
            endif
         enddo
      else                      !@@@@@@@@@@@@@@@@@@@@
         nseq2_=L_fr
         aL=min(nseq1,nseq2_)
         idel=aL/2.5            !minimum size of considered fragment
         if(idel.le.fra_min1)idel=fra_min1
         n1=-nseq2_+idel
         n2=nseq1-idel
         GL_max=0
         do ishift=n1,n2
            L=0
            do j=1,nseq2
               invmap(j)=-1
            enddo
            do j=1,nseq2_
               i=j+ishift
               if(i.ge.1.and.i.le.nseq1)then
                  L=L+1
                  invmap(ifr(j))=i
               endif
            enddo
            if(L.ge.idel)then
               call get_GL(GL)
               if(GL.gt.GL_max)then
                  GL_max=GL
                  do i=1,nseq2
                     invmap_i(i)=invmap(i)
                  enddo
               endif
            endif
         enddo
      endif
      
      return
      end

**************************************************************
*     smooth the secondary structure assignment
**************************************************************
      subroutine smooth
      PARAMETER(nmax=5000)
      common/sec/isec(nmax),jsec(nmax)
      common/length/nseq1,nseq2

***   smooth single -------------->
***   --x-- => -----
      do i=1,nseq1
         if(isec(i).eq.2.or.isec(i).eq.4)then
            j=isec(i)
            if(isec(i-2).ne.j)then
               if(isec(i-1).ne.j)then
                  if(isec(i+1).ne.j)then
                     if(isec(i+1).ne.j)then
                        isec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).eq.2.or.jsec(i).eq.4)then
            j=jsec(i)
            if(jsec(i-2).ne.j)then
               if(jsec(i-1).ne.j)then
                  if(jsec(i+1).ne.j)then
                     if(jsec(i+1).ne.j)then
                        jsec(i)=1
                     endif
                  endif
               endif
            endif
         endif
      enddo

***   smooth double -------------->
***   --xx-- => ------
      do i=1,nseq1
         if(isec(i).ne.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
         if(isec(i+3).eq.2)then
         if(isec(i+4).ne.2)then
         if(isec(i+5).ne.2)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(isec(i).ne.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
         if(isec(i+3).eq.4)then
         if(isec(i+4).ne.4)then
         if(isec(i+5).ne.4)then
            isec(i+2)=1
            isec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).ne.2)then
         if(jsec(i+1).ne.2)then
         if(jsec(i+2).eq.2)then
         if(jsec(i+3).eq.2)then
         if(jsec(i+4).ne.2)then
         if(jsec(i+5).ne.2)then
            jsec(i+2)=1
            jsec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif

         if(jsec(i).ne.4)then
         if(jsec(i+1).ne.4)then
         if(jsec(i+2).eq.4)then
         if(jsec(i+3).eq.4)then
         if(jsec(i+4).ne.4)then
         if(jsec(i+5).ne.4)then
            jsec(i+2)=1
            jsec(i+3)=1
         endif
         endif
         endif
         endif
         endif
         endif
      enddo

***   connect -------------->
***   x-x => xxx
      do i=1,nseq1
         if(isec(i).eq.2)then
         if(isec(i+1).ne.2)then
         if(isec(i+2).eq.2)then
            isec(i+1)=2
         endif
         endif
         endif

         if(isec(i).eq.4)then
         if(isec(i+1).ne.4)then
         if(isec(i+2).eq.4)then
            isec(i+1)=4
         endif
         endif
         endif
      enddo
      do i=1,nseq2
         if(jsec(i).eq.2)then
         if(jsec(i+1).ne.2)then
         if(jsec(i+2).eq.2)then
            jsec(i+1)=2
         endif
         endif
         endif

         if(jsec(i).eq.4)then
         if(jsec(i+1).ne.4)then
         if(jsec(i+2).eq.4)then
            jsec(i+1)=4
         endif
         endif
         endif
      enddo

      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function diszy(i,i1,i2)
      PARAMETER(nmax=5000)
      COMMON/BACKBONE/XA(3,nmax,0:1)
      diszy=sqrt((xa(1,i1,i)-xa(1,i2,i))**2
     &     +(xa(2,i1,i)-xa(2,i2,i))**2
     &     +(xa(3,i1,i)-xa(3,i2,i))**2)
      return
      end

*************************************************************
*     assign secondary structure:
*************************************************************
      function make_sec(dis13,dis14,dis15,dis24,dis25,dis35)
      make_sec=1
      delta=2.1
      if(abs(dis15-6.37).lt.delta)then
         if(abs(dis14-5.18).lt.delta)then
            if(abs(dis25-5.18).lt.delta)then
               if(abs(dis13-5.45).lt.delta)then
                  if(abs(dis24-5.45).lt.delta)then
                     if(abs(dis35-5.45).lt.delta)then
                        make_sec=2 !helix
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      delta=1.42
      if(abs(dis15-13).lt.delta)then
         if(abs(dis14-10.4).lt.delta)then
            if(abs(dis25-10.4).lt.delta)then
               if(abs(dis13-6.1).lt.delta)then
                  if(abs(dis24-6.1).lt.delta)then
                     if(abs(dis35-6.1).lt.delta)then
                        make_sec=4 !strand
                        return
                     endif
                  endif
               endif
            endif
         endif
      endif
      if(dis15.lt.8)then
         make_sec=3
      endif

      return
      end

****************************************************************
*     quickly calculate TM-score with given invmap(i) in 3 iterations
****************************************************************
      subroutine get_GL(GL)
      PARAMETER(nmax=5000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)
      common/d00/d00,d002

      dimension xo1(nmax),yo1(nmax),zo1(nmax)
      dimension xo2(nmax),yo2(nmax),zo2(nmax)
      dimension dis2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
            xo1(n_al)=xa(1,i,0)
            yo1(n_al)=xa(2,i,0)
            zo1(n_al)=xa(3,i,0)
            xo2(n_al)=xa(1,j,1)
            yo2(n_al)=xa(2,j,1)
            zo2(n_al)=xa(3,j,1)
         endif
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      GL=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         GL=GL+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration------------->
      d002t=d002
 21   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3.and.n_al.gt.3)then
         d002t=d002t+.5
         goto 21
      endif
      L=j
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
      G2=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G2=G2+1/(1+dis2(i)/(d0**2))
      enddo
ccc   for next iteration------------->
      d002t=d002+1
 22   j=0
      do i=1,n_al
         if(dis2(i).le.d002t)then
            j=j+1
            r_1(1,j)=xo1(i)
            r_1(2,j)=yo1(i)
            r_1(3,j)=zo1(i)
            r_2(1,j)=xo2(i)
            r_2(2,j)=yo2(i)
            r_2(3,j)=zo2(i)
         endif
      enddo
      if(j.lt.3.and.n_al.gt.3)then
         d002t=d002t+.5
         goto 22
      endif
      L=j
      call u3b(w,r_1,r_2,L,1,rms,u,t,ier) !u rotate r_1 to r_2
      G3=0
      do i=1,n_al
         xx=t(1)+u(1,1)*xo1(i)+u(1,2)*yo1(i)+u(1,3)*zo1(i)
         yy=t(2)+u(2,1)*xo1(i)+u(2,2)*yo1(i)+u(2,3)*zo1(i)
         zz=t(3)+u(3,1)*xo1(i)+u(3,2)*yo1(i)+u(3,3)*zo1(i)
         dis2(i)=(xx-xo2(i))**2+(yy-yo2(i))**2+(zz-zo2(i))**2
         G3=G3+1/(1+dis2(i)/(d0**2))
      enddo
      if(G2.gt.GL)GL=G2
      if(G3.gt.GL)GL=G3

c^^^^^^^^^^^^^^^^ GL done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

****************************************************************
*     with invmap(i) calculate TM-score and martix score(i,j) for rotation 
****************************************************************
      subroutine get_score
      PARAMETER(nmax=5000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for TM-score:
            xtm1(n_al)=xa(1,i,0) !for TM-score
            ytm1(n_al)=xa(2,i,0)
            ztm1(n_al)=xa(3,i,0)
            xtm2(n_al)=xa(1,j,1)
            ytm2(n_al)=xa(2,j,1)
            ztm2(n_al)=xa(3,j,1)
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
         endif
      enddo
***   calculate TM-score for the given alignment----------->
      d0_input=d0
      call TMscore8_search(d0_input,n_al,xtm1,ytm1,ztm1,n1,
     &     n_al,xtm2,ytm2,ztm2,n2,TM,Rcomm,Lcomm) !simplified search engine
      TM=TM*n_al/anseq          !TM-score
***   calculate score matrix score(i,j)------------------>
      do i=1,n_al
         r_2(1,i)=xtm1(i)
         r_2(2,i)=ytm1(i)
         r_2(3,i)=ztm1(i)
      enddo
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d0**2)
         enddo
      enddo
      
c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end

****************************************************************
*     with invmap(i) calculate score(i,j) using RMSD rotation 
****************************************************************
      subroutine get_score1
      PARAMETER(nmax=5000)
      common/length/nseq1,nseq2
      COMMON/BACKBONE/XA(3,nmax,0:1)
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      common/zscore/zrms,n_al,rmsd_al
      common/d0/d0,anseq
      common/d0min/d0_min
      dimension xtm1(nmax),ytm1(nmax),ztm1(nmax)
      dimension xtm2(nmax),ytm2(nmax),ztm2(nmax)
      common/TM/TM,TMmax
      common/n1n2/n1(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

c     calculate RMSD between aligned structures and rotate the structures -->
      n_al=0
      do j=1,NSEQ2
         i=invmap(j)            !j aligned to i
         if(i.gt.0)then
            n_al=n_al+1
ccc   for rotation matrix:
            r_1(1,n_al)=xa(1,i,0)
            r_1(2,n_al)=xa(2,i,0)
            r_1(3,n_al)=xa(3,i,0)
            r_2(1,n_al)=xa(1,j,1)
            r_2(2,n_al)=xa(2,j,1)
            r_2(3,n_al)=xa(3,j,1)
         endif
      enddo
***   calculate score matrix score(i,j)------------------>
      call u3b(w,r_1,r_2,n_al,1,rms,u,t,ier) !u rotate r_1 to r_2
      d01=d0+1.5
      if(d01.lt.d0_min)d01=d0_min
      d02=d01*d01
      do i=1,nseq1
         xx=t(1)+u(1,1)*xa(1,i,0)+u(1,2)*xa(2,i,0)+u(1,3)*xa(3,i,0)
         yy=t(2)+u(2,1)*xa(1,i,0)+u(2,2)*xa(2,i,0)+u(2,3)*xa(3,i,0)
         zz=t(3)+u(3,1)*xa(1,i,0)+u(3,2)*xa(2,i,0)+u(3,3)*xa(3,i,0)
         do j=1,nseq2
            dd=(xx-xa(1,j,1))**2+(yy-xa(2,j,1))**2+(zz-xa(3,j,1))**2
            score(i,j)=1/(1+dd/d02)
         enddo
      enddo

c^^^^^^^^^^^^^^^^ score(i,j) done ^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      end


*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
*** dis<8, simplified search engine
      subroutine TMscore8_search(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)
      dimension iL0(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        k=0
        do i=1,iL_max,40        !this is the simplification!
          k=k+1
          iL0(k)=i
        enddo
        if(iL0(k).lt.iL_max)then
          k=k+1
          iL0(k)=iL_max
        endif
        n_shift=k
        do 300 i_shift=1,n_shift
           iL=iL0(i_shift)
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun8       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun8    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

******** return the final rotation ****************
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
***   dis<8, but same search engine
      subroutine TMscore8(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max    !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun8       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun8    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

******** return the final rotation ****************
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun8
      PARAMETER(nmax=5000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score,score_max
      common/d8/d8

      d_tmp=d
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d_tmp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         if(dis.le.d8)then
            score_sum=score_sum+1/(1+(dis/d0)**2)
         endif
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         goto 21
      endif
      score=score_sum/float(nseqB) !TM-score

      return
      end

*************************************************************************
*************************************************************************
*     This is a subroutine to compare two structures and find the 
*     superposition that has the maximum TM-score.
*
*     L1--Length of the first structure
*     (x1(i),y1(i),z1(i))--coordinates of i'th residue at the first structure
*     n1(i)--Residue sequence number of i'th residue at the first structure
*     L2--Length of the second structure
*     (x2(i),y2(i),z2(i))--coordinates of i'th residue at the second structure
*     n2(i)--Residue sequence number of i'th residue at the second structure
*     TM--TM-score of the comparison
*     Rcomm--RMSD of two structures in the common aligned residues
*     Lcomm--Length of the common aligned regions
*
*     Note: 
*     1, Always put native as the second structure, by which TM-score
*        is normalized.
*     2, The returned (x1(i),y1(i),z1(i)) are the rotated structure after
*        TM-score superposition.
*************************************************************************
*************************************************************************
***  normal TM-score:
      subroutine TMscore(dx,L1,x1,y1,z1,n1,L2,x2,y2,z2,n2,
     &     TM,Rcomm,Lcomm)
      PARAMETER(nmax=5000)
      common/stru/xt(nmax),yt(nmax),zt(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/d0min/d0_min
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      dimension k_ali(nmax),k_ali0(nmax)
      dimension L_ini(100),iq(nmax)
      common/scores/score
      double precision score,score_max
      dimension xa(nmax),ya(nmax),za(nmax)

      dimension x1(nmax),y1(nmax),z1(nmax),n1(nmax)
      dimension x2(nmax),y2(nmax),z2(nmax),n2(nmax)

ccc   RMSD:
      double precision r_1(3,nmax),r_2(3,nmax),r_3(3,nmax),w(nmax)
      double precision u(3,3),t(3),rms,drms !armsd is real
      data w /nmax*1.0/
ccc   

********* convert input data ***************************
*     because L1=L2 in this special case---------->
      nseqA=L1
      nseqB=L2
      do i=1,nseqA
         xa(i)=x1(i)
         ya(i)=y1(i)
         za(i)=z1(i)
         nresA(i)=n1(i)
         xb(i)=x2(i)
         yb(i)=y2(i)
         zb(i)=z2(i)
         nresB(i)=n2(i)
         iA(i)=i
         iB(i)=i
      enddo
      n_ali=L1                  !number of aligned residues
      Lcomm=L1

************/////
*     parameters:
*****************
***   d0------------->
c      d0=1.24*(nseqB-15)**(1.0/3.0)-1.8
      d0=dx
      if(d0.lt.d0_min)d0=d0_min
***   d0_search ----->
      d0_search=d0
      if(d0_search.gt.8)d0_search=8
      if(d0_search.lt.4.5)d0_search=4.5
***   iterative parameters ----->
      n_it=20                   !maximum number of iterations
      d_output=5                !for output alignment
      n_init_max=6              !maximum number of L_init
      n_init=0
      L_ini_min=4
      if(n_ali.lt.4)L_ini_min=n_ali
      do i=1,n_init_max-1
         n_init=n_init+1
         L_ini(n_init)=n_ali/2**(n_init-1)
         if(L_ini(n_init).le.L_ini_min)then
            L_ini(n_init)=L_ini_min
            goto 402
         endif
      enddo
      n_init=n_init+1
      L_ini(n_init)=L_ini_min
 402  continue

******************************************************************
*     find the maximum score starting from local structures superposition
******************************************************************
      score_max=-1              !TM-score
      do 333 i_init=1,n_init
        L_init=L_ini(i_init)
        iL_max=n_ali-L_init+1
        do 300 iL=1,iL_max      !on aligned residues, [1,nseqA]
           LL=0
           ka=0
           do i=1,L_init
              k=iL+i-1          ![1,n_ali] common aligned
              r_1(1,i)=xa(iA(k))
              r_1(2,i)=ya(iA(k))
              r_1(3,i)=za(iA(k))
              r_2(1,i)=xb(iB(k))
              r_2(2,i)=yb(iB(k))
              r_2(3,i)=zb(iB(k))
              LL=LL+1
              ka=ka+1
              k_ali(ka)=k
           enddo
           call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
           if(i_init.eq.1)then  !global superposition
              armsd=dsqrt(rms/LL)
              Rcomm=armsd
           endif
           do j=1,nseqA
              xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
              yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
              zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
           enddo
           d=d0_search-1
           call score_fun       !init, get scores, n_cut+i_ali(i) for iteration
           if(score_max.lt.score)then
              score_max=score
              ka0=ka
              do i=1,ka0
                 k_ali0(i)=k_ali(i)
              enddo
           endif
***   iteration for extending ---------------------------------->
           d=d0_search+1
           do 301 it=1,n_it
              LL=0
              ka=0
              do i=1,n_cut
                 m=i_ali(i)     ![1,n_ali]
                 r_1(1,i)=xa(iA(m))
                 r_1(2,i)=ya(iA(m))
                 r_1(3,i)=za(iA(m))
                 r_2(1,i)=xb(iB(m))
                 r_2(2,i)=yb(iB(m))
                 r_2(3,i)=zb(iB(m))
                 ka=ka+1
                 k_ali(ka)=m
                 LL=LL+1
              enddo
              call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
              do j=1,nseqA
                 xt(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
                 yt(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
                 zt(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
              enddo
              call score_fun    !get scores, n_cut+i_ali(i) for iteration
              if(score_max.lt.score)then
                 score_max=score
                 ka0=ka
                 do i=1,ka
                    k_ali0(i)=k_ali(i)
                 enddo
              endif
              if(it.eq.n_it)goto 302
              if(n_cut.eq.ka)then
                 neq=0
                 do i=1,n_cut
                    if(i_ali(i).eq.k_ali(i))neq=neq+1
                 enddo
                 if(n_cut.eq.neq)goto 302
              endif
 301       continue             !for iteration
 302       continue
 300    continue                !for shift
 333  continue                  !for initial length, L_ali/M

******** return the final rotation ****************
      LL=0
      do i=1,ka0
         m=k_ali0(i)            !record of the best alignment
         r_1(1,i)=xa(iA(m))
         r_1(2,i)=ya(iA(m))
         r_1(3,i)=za(iA(m))
         r_2(1,i)=xb(iB(m))
         r_2(2,i)=yb(iB(m))
         r_2(3,i)=zb(iB(m))
         LL=LL+1
      enddo
      call u3b(w,r_1,r_2,LL,1,rms,u,t,ier) !u rotate r_1 to r_2
      do j=1,nseqA
         x1(j)=t(1)+u(1,1)*xa(j)+u(1,2)*ya(j)+u(1,3)*za(j)
         y1(j)=t(2)+u(2,1)*xa(j)+u(2,2)*ya(j)+u(2,3)*za(j)
         z1(j)=t(3)+u(3,1)*xa(j)+u(3,2)*ya(j)+u(3,3)*za(j)
      enddo
      TM=score_max

c^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
      return
      END

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1, collect those residues with dis<d;
c     2, calculate score_GDT, score_maxsub, score_TM
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine score_fun
      PARAMETER(nmax=5000)

      common/stru/xa(nmax),ya(nmax),za(nmax),xb(nmax),yb(nmax),zb(nmax)
      common/nres/nresA(nmax),nresB(nmax),nseqA,nseqB
      common/para/d,d0
      common/align/n_ali,iA(nmax),iB(nmax)
      common/nscore/i_ali(nmax),n_cut ![1,n_ali],align residues for the score
      common/scores/score
      double precision score

      d_tmp=d
 21   n_cut=0                   !number of residue-pairs dis<d, for iteration
      score_sum=0               !TMscore
      do k=1,n_ali
         i=iA(k)                ![1,nseqA] reoder number of structureA
         j=iB(k)                ![1,nseqB]
         dis=sqrt((xa(i)-xb(j))**2+(ya(i)-yb(j))**2+(za(i)-zb(j))**2)
         if(dis.lt.d_tmp)then
            n_cut=n_cut+1
            i_ali(n_cut)=k      ![1,n_ali], mark the residue-pairs in dis<d
         endif
         score_sum=score_sum+1/(1+(dis/d0)**2)
      enddo
      if(n_cut.lt.3.and.n_ali.gt.3)then
         d_tmp=d_tmp+.5
         goto 21
      endif
      score=score_sum/float(nseqB) !TM-score

      return
      end

********************************************************************
*     Dynamic programming for alignment.
*     Input: score(i,j), and gap_open
*     Output: invmap(j)
*     
*     Please note this subroutine is not a correct implementation of 
*     the N-W dynamic programming because the score tracks back only 
*     one layer of the matrix. This code was exploited in TM-align 
*     because it is about 1.5 times faster than a complete N-W code
*     and does not influence much the final structure alignment result.
********************************************************************
      SUBROUTINE DP(NSEQ1,NSEQ2)
      PARAMETER(nmax=5000)
      LOGICAL*1 DIR
      common/dpc/score(nmax,nmax),gap_open,invmap(nmax)
      dimension DIR(0:nmax,0:nmax),VAL(0:nmax,0:nmax)
      REAL H,V
      
***   initialize the matrix:
      val(0,0)=0
      do i=1,nseq1
        dir(i,0)=.false.
        val(i,0)=0
      enddo
      do j=1,nseq2
        dir(0,j)=.false.
        val(0,j)=0
        invmap(j)=-1
      enddo

***   decide matrix and path:
      DO j=1,NSEQ2
        DO i=1,NSEQ1
          D=VAL(i-1,j-1)+SCORE(i,j)
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          
          IF((D.GE.H).AND.(D.GE.V)) THEN
            DIR(I,J)=.true.
            VAL(i,j)=D
          ELSE
            DIR(I,J)=.false.
            if(V.GE.H)then
              val(i,j)=v
            else
              val(i,j)=h
            end if
          ENDIF
        ENDDO
      ENDDO
      
***   extract the alignment:
      i=NSEQ1
      j=NSEQ2
      DO WHILE((i.GT.0).AND.(j.GT.0))
        IF(DIR(i,j))THEN
          invmap(j)=i
          i=i-1
          j=j-1
        ELSE
          H=VAL(i-1,j)
          if(DIR(i-1,j))H=H+GAP_OPEN
          V=VAL(i,j-1)
          if(DIR(i,j-1))V=V+GAP_OPEN
          IF(V.GE.H) THEN
            j=j-1
          ELSE
            i=i-1
          ENDIF
        ENDIF
      ENDDO
      
c^^^^^^^^^^^^^^^Dynamical programming done ^^^^^^^^^^^^^^^^^^^
      return
      END

cccccccccccccccc Calculate sum of (r_d-r_m)^2 cccccccccccccccccccccccccc
c  w    - w(m) is weight for atom pair  c m           (given)
c  x    - x(i,m) are coordinates of atom c m in set x       (given)
c  y    - y(i,m) are coordinates of atom c m in set y       (given)
c  n    - n is number of atom pairs                         (given)
c  mode  - 0:calculate rms only                             (given)
c          1:calculate rms,u,t                              (takes longer)
c  rms   - sum of w*(ux+t-y)**2 over all atom pairs         (result)
c  u    - u(i,j) is   rotation  matrix for best superposition  (result)
c  t    - t(i)   is translation vector for best superposition  (result)
c  ier  - 0: a unique optimal superposition has been determined(result)
c       -1: superposition is not unique but optimal
c       -2: no result obtained because of negative weights w
c           or all weights equal to zero.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine u3b(w, x, y, n, mode, rms, u, t, ier)
      integer ip(9), ip2312(4), i, j, k, l, m1, m, ier, n, mode
      double precision w(1), x(3, 1), y(3, 1), u(3, 3), t(3), rms, sigma
      double precision r(3, 3), xc(3), yc(3), wc, a(3, 3), b(3, 3), e0, 
     &e(3), e1, e2, e3, d, spur, det, cof, h, g, cth, sth, sqrth, p, tol
     &, rr(6), rr1, rr2, rr3, rr4, rr5, rr6, ss(6), ss1, ss2, ss3, ss4, 
     &ss5, ss6, zero, one, two, three, sqrt3
      equivalence (rr(1), rr1), (rr(2), rr2), (rr(3), rr3), (rr(4), rr4)
     &, (rr(5), rr5), (rr(6), rr6), (ss(1), ss1), (ss(2), ss2), (ss(3), 
     &ss3), (ss(4), ss4), (ss(5), ss5), (ss(6), ss6), (e(1), e1), (e(2)
     &, e2), (e(3), e3)
      data sqrt3 / 1.73205080756888d+00 /
      data tol / 1.0d-2 /
      data zero / 0.0d+00 /
      data one / 1.0d+00 /
      data two / 2.0d+00 /
      data three / 3.0d+00 /
      data ip / 1, 2, 4, 2, 3, 5, 4, 5, 6 /
      data ip2312 / 2, 3, 1, 2 /
c 156 "rms.for"
      wc = zero
      rms = 0.0
      e0 = zero
      do 1 i = 1, 3
      xc(i) = zero
      yc(i) = zero
      t(i) = 0.0
      do 1 j = 1, 3
      d = zero
      if (i .eq. j) d = one
      u(i,j) = d
      a(i,j) = d
    1 r(i,j) = zero
      ier = -1
c**** DETERMINE CENTROIDS OF BOTH VECTOR SETS X AND Y
c 170 "rms.for"
      if (n .lt. 1) return 
c 172 "rms.for"
      ier = -2
      do 2 m = 1, n
      if (w(m) .lt. 0.0) return 
      wc = wc + w(m)
      do 2 i = 1, 3
      xc(i) = xc(i) + (w(m) * x(i,m))
    2 yc(i) = yc(i) + (w(m) * y(i,m))
      if (wc .le. zero) return 
      do 3 i = 1, 3
      xc(i) = xc(i) / wc
c**** DETERMINE CORRELATION MATRIX R BETWEEN VECTOR SETS Y AND X
c 182 "rms.for"
    3 yc(i) = yc(i) / wc
c 184 "rms.for"
      do 4 m = 1, n
      do 4 i = 1, 3
      e0 = e0 + (w(m) * (((x(i,m) - xc(i)) ** 2) + ((y(i,m) - yc(i)) ** 
     &2)))
c 187 "rms.for"
      d = w(m) * (y(i,m) - yc(i))
      do 4 j = 1, 3
c**** CALCULATE DETERMINANT OF R(I,J)
c 189 "rms.for"
    4 r(i,j) = r(i,j) + (d * (x(j,m) - xc(j)))
c 191 "rms.for"
      det = ((r(1,1) * ((r(2,2) * r(3,3)) - (r(2,3) * r(3,2)))) - (r(1,2
     &) * ((r(2,1) * r(3,3)) - (r(2,3) * r(3,1))))) + (r(1,3) * ((r(2,1)
     & * r(3,2)) - (r(2,2) * r(3,1))))
c**** FORM UPPER TRIANGLE OF TRANSPOSED(R)*R
c 194 "rms.for"
      sigma = det
c 196 "rms.for"
      m = 0
      do 5 j = 1, 3
      do 5 i = 1, j
      m = m + 1
c***************** EIGENVALUES *****************************************
c**** FORM CHARACTERISTIC CUBIC  X**3-3*SPUR*X**2+3*COF*X-DET=0
c 200 "rms.for"
    5 rr(m) = ((r(1,i) * r(1,j)) + (r(2,i) * r(2,j))) + (r(3,i) * r(3,j)
     &)
c 203 "rms.for"
      spur = ((rr1 + rr3) + rr6) / three
      cof = ((((((rr3 * rr6) - (rr5 * rr5)) + (rr1 * rr6)) - (rr4 * rr4)
     &) + (rr1 * rr3)) - (rr2 * rr2)) / three
c 205 "rms.for"
      det = det * det
      do 6 i = 1, 3
    6 e(i) = spur
c**** REDUCE CUBIC TO STANDARD FORM Y**3-3HY+2G=0 BY PUTTING X=Y+SPUR
c 208 "rms.for"
      if (spur .le. zero) goto 40
c 210 "rms.for"
      d = spur * spur
      h = d - cof
c**** SOLVE CUBIC. ROOTS ARE E1,E2,E3 IN DECREASING ORDER
c 212 "rms.for"
      g = (((spur * cof) - det) / two) - (spur * h)
c 214 "rms.for"
      if (h .le. zero) goto 8
      sqrth = dsqrt(h)
      d = ((h * h) * h) - (g * g)
      if (d .lt. zero) d = zero
      d = datan2(dsqrt(d),- g) / three
      cth = sqrth * dcos(d)
      sth = (sqrth * sqrt3) * dsin(d)
      e1 = (spur + cth) + cth
      e2 = (spur - cth) + sth
      e3 = (spur - cth) - sth
c.....HANDLE SPECIAL CASE OF 3 IDENTICAL ROOTS
c 224 "rms.for"
      if (mode .eq. 0) goto 50
      if (mode .ne. 0) goto 10
c**************** EIGENVECTORS *****************************************
c 226 "rms.for"
    8 if (mode .eq. 0) goto 50
      if (mode .ne. 0) goto 30
c 228 "rms.for"
   10 do 15 l = 1, 3, 2
      d = e(l)
      ss1 = ((d - rr3) * (d - rr6)) - (rr5 * rr5)
      ss2 = ((d - rr6) * rr2) + (rr4 * rr5)
      ss3 = ((d - rr1) * (d - rr6)) - (rr4 * rr4)
      ss4 = ((d - rr3) * rr4) + (rr2 * rr5)
      ss5 = ((d - rr1) * rr5) + (rr2 * rr4)
      ss6 = ((d - rr1) * (d - rr3)) - (rr2 * rr2)
      j = 1
      if (dabs(ss1) .ge. dabs(ss3)) goto 12
      j = 2
      if (dabs(ss3) .ge. dabs(ss6)) goto 13
   11 j = 3
      goto 13
   12 if (dabs(ss1) .lt. dabs(ss6)) goto 11
   13 d = zero
      j = 3 * (j - 1)
      do 14 i = 1, 3
      k = ip(i + j)
      a(i,l) = ss(k)
   14 d = d + (ss(k) * ss(k))
      if (d .gt. zero) d = one / dsqrt(d)
      do 15 i = 1, 3
   15 a(i,l) = a(i,l) * d
      d = ((a(1,1) * a(1,3)) + (a(2,1) * a(2,3))) + (a(3,1) * a(3,3))
      m1 = 3
      m = 1
      if ((e1 - e2) .gt. (e2 - e3)) goto 16
      m1 = 1
      m = 3
   16 p = zero
      do 17 i = 1, 3
      a(i,m1) = a(i,m1) - (d * a(i,m))
   17 p = p + (a(i,m1) ** 2)
      if (p .le. tol) goto 19
      p = one / dsqrt(p)
      do 18 i = 1, 3
   18 a(i,m1) = a(i,m1) * p
      goto 21
   19 p = one
      do 20 i = 1, 3
      if (p .lt. dabs(a(i,m))) goto 20
      p = dabs(a(i,m))
      j = i
   20 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((a(k,m) ** 2) + (a(l,m) ** 2))
      if (p .le. tol) goto 40
      a(j,m1) = zero
      a(k,m1) = - (a(l,m) / p)
      a(l,m1) = a(k,m) / p
   21 a(1,2) = (a(2,3) * a(3,1)) - (a(2,1) * a(3,3))
      a(2,2) = (a(3,3) * a(1,1)) - (a(3,1) * a(1,3))
c****************** ROTATION MATRIX ************************************
c 282 "rms.for"
      a(3,2) = (a(1,3) * a(2,1)) - (a(1,1) * a(2,3))
c 284 "rms.for"
   30 do 32 l = 1, 2
      d = zero
      do 31 i = 1, 3
      b(i,l) = ((r(i,1) * a(1,l)) + (r(i,2) * a(2,l))) + (r(i,3) * a(3,l
     &))
c 288 "rms.for"
   31 d = d + (b(i,l) ** 2)
      if (d .gt. zero) d = one / dsqrt(d)
      do 32 i = 1, 3
   32 b(i,l) = b(i,l) * d
      d = ((b(1,1) * b(1,2)) + (b(2,1) * b(2,2))) + (b(3,1) * b(3,2))
      p = zero
      do 33 i = 1, 3
      b(i,2) = b(i,2) - (d * b(i,1))
   33 p = p + (b(i,2) ** 2)
      if (p .le. tol) goto 35
      p = one / dsqrt(p)
      do 34 i = 1, 3
   34 b(i,2) = b(i,2) * p
      goto 37
   35 p = one
      do 36 i = 1, 3
      if (p .lt. dabs(b(i,1))) goto 36
      p = dabs(b(i,1))
      j = i
   36 continue
      k = ip2312(j)
      l = ip2312(j + 1)
      p = dsqrt((b(k,1) ** 2) + (b(l,1) ** 2))
      if (p .le. tol) goto 40
      b(j,2) = zero
      b(k,2) = - (b(l,1) / p)
      b(l,2) = b(k,1) / p
   37 b(1,3) = (b(2,1) * b(3,2)) - (b(2,2) * b(3,1))
      b(2,3) = (b(3,1) * b(1,2)) - (b(3,2) * b(1,1))
      b(3,3) = (b(1,1) * b(2,2)) - (b(1,2) * b(2,1))
      do 39 i = 1, 3
      do 39 j = 1, 3
c****************** TRANSLATION VECTOR *********************************
c 320 "rms.for"
   39 u(i,j) = ((b(i,1) * a(j,1)) + (b(i,2) * a(j,2))) + (b(i,3) * a(j,3
     &))
   40 do 41 i = 1, 3
c********************** RMS ERROR **************************************
c 323 "rms.for"
   41 t(i) = ((yc(i) - (u(i,1) * xc(1))) - (u(i,2) * xc(2))) - (u(i,3)
     & * xc(3))
   50 do 51 i = 1, 3
      if (e(i) .lt. zero) e(i) = zero
   51 e(i) = dsqrt(e(i))
      ier = 0
      if (e2 .le. (e1 * 1.0d-05)) ier = -1
      d = e3
      if (sigma .ge. 0.0) goto 52
      d = - d
      if ((e2 - e3) .le. (e1 * 1.0d-05)) ier = -1
   52 d = (d + e2) + e1
      rms = (e0 - d) - d
      if (rms .lt. 0.0) rms = 0.0
      return 
c.....END U3B...........................................................
c----------------------------------------------------------
c                       THE END
c----------------------------------------------------------
c 338 "rms.for"
      end

