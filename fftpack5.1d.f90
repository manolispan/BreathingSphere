subroutine c1f2kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F2KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    13 May 2013
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,2)
  real ( kind = 8 ) ch(in2,l1,2,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,1,2)

  if ( ido <= 1 ) then

    sn = 1.0D+00 / real ( 2 * l1, kind = 8 )

    if ( na == 1 ) then

      do k = 1, l1
        ch(1,k,1,1) = sn*(cc(1,k,1,1)+cc(1,k,1,2))
        ch(1,k,2,1) = sn*(cc(1,k,1,1)-cc(1,k,1,2))
        ch(2,k,1,1) = sn*(cc(2,k,1,1)+cc(2,k,1,2))
        ch(2,k,2,1) = sn*(cc(2,k,1,1)-cc(2,k,1,2))
      end do

    else

      do k = 1, l1
        chold1 = sn*(cc(1,k,1,1)+cc(1,k,1,2))
        cc(1,k,1,2) = sn*(cc(1,k,1,1)-cc(1,k,1,2))
        cc(1,k,1,1) = chold1
        chold2 = sn*(cc(2,k,1,1)+cc(2,k,1,2))
        cc(2,k,1,2) = sn*(cc(2,k,1,1)-cc(2,k,1,2))
        cc(2,k,1,1) = chold2
      end do

    end if

  else

    do k = 1, l1
      ch(1,k,1,1) = cc(1,k,1,1)+cc(1,k,1,2)
      ch(1,k,2,1) = cc(1,k,1,1)-cc(1,k,1,2)
      ch(2,k,1,1) = cc(2,k,1,1)+cc(2,k,1,2)
      ch(2,k,2,1) = cc(2,k,1,1)-cc(2,k,1,2)
    end do

    do i = 2, ido
      do k = 1, l1
        ch(1,k,1,i) = cc(1,k,i,1)+cc(1,k,i,2)
        tr2 = cc(1,k,i,1)-cc(1,k,i,2)
        ch(2,k,1,i) = cc(2,k,i,1)+cc(2,k,i,2)
        ti2 = cc(2,k,i,1)-cc(2,k,i,2)
        ch(2,k,2,i) = wa(i,1,1)*ti2-wa(i,1,2)*tr2
        ch(1,k,2,i) = wa(i,1,1)*tr2+wa(i,1,2)*ti2
      end do
    end do

  end if

  return
end


subroutine c1f3kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F3KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    13 May 2013
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,3)
  real ( kind = 8 ) ch(in2,l1,3,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ), parameter :: taui = -0.866025403784439D+00
  real ( kind = 8 ), parameter :: taur = -0.5D+00
  real ( kind = 8 ) ti2
  real ( kind = 8 ) tr2
  real ( kind = 8 ) wa(ido,2,2)

  if ( ido <= 1 ) then

    sn = 1.0D+00 / real ( 3 * l1, kind = 8 )

    if ( na /= 1 ) then

      do k = 1, l1
        tr2 = cc(1,k,1,2)+cc(1,k,1,3)
        cr2 = cc(1,k,1,1)+taur*tr2
        cc(1,k,1,1) = sn*(cc(1,k,1,1)+tr2)
        ti2 = cc(2,k,1,2)+cc(2,k,1,3)
        ci2 = cc(2,k,1,1)+taur*ti2
        cc(2,k,1,1) = sn*(cc(2,k,1,1)+ti2)
        cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
        ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
        cc(1,k,1,2) = sn*(cr2-ci3)
        cc(1,k,1,3) = sn*(cr2+ci3)
        cc(2,k,1,2) = sn*(ci2+cr3)
        cc(2,k,1,3) = sn*(ci2-cr3)
      end do

    else

      do k = 1, l1
        tr2 = cc(1,k,1,2)+cc(1,k,1,3)
        cr2 = cc(1,k,1,1)+taur*tr2
        ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2)
        ti2 = cc(2,k,1,2)+cc(2,k,1,3)
        ci2 = cc(2,k,1,1)+taur*ti2
        ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2)
        cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
        ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
        ch(1,k,2,1) = sn*(cr2-ci3)
        ch(1,k,3,1) = sn*(cr2+ci3)
        ch(2,k,2,1) = sn*(ci2+cr3)
        ch(2,k,3,1) = sn*(ci2-cr3)
      end do

    end if

  else

    do k = 1, l1
      tr2 = cc(1,k,1,2)+cc(1,k,1,3)
      cr2 = cc(1,k,1,1)+taur*tr2
      ch(1,k,1,1) = cc(1,k,1,1)+tr2
      ti2 = cc(2,k,1,2)+cc(2,k,1,3)
      ci2 = cc(2,k,1,1)+taur*ti2
      ch(2,k,1,1) = cc(2,k,1,1)+ti2
      cr3 = taui*(cc(1,k,1,2)-cc(1,k,1,3))
      ci3 = taui*(cc(2,k,1,2)-cc(2,k,1,3))
      ch(1,k,2,1) = cr2-ci3
      ch(1,k,3,1) = cr2+ci3
      ch(2,k,2,1) = ci2+cr3
      ch(2,k,3,1) = ci2-cr3
    end do

    do i = 2, ido
      do k = 1, l1
        tr2 = cc(1,k,i,2)+cc(1,k,i,3)
        cr2 = cc(1,k,i,1)+taur*tr2
        ch(1,k,1,i) = cc(1,k,i,1)+tr2
        ti2 = cc(2,k,i,2)+cc(2,k,i,3)
        ci2 = cc(2,k,i,1)+taur*ti2
        ch(2,k,1,i) = cc(2,k,i,1)+ti2
        cr3 = taui*(cc(1,k,i,2)-cc(1,k,i,3))
        ci3 = taui*(cc(2,k,i,2)-cc(2,k,i,3))
        dr2 = cr2-ci3
        dr3 = cr2+ci3
        di2 = ci2+cr3
        di3 = ci2-cr3
        ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
        ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
        ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
        ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
      end do
    end do

  end if

  return
end


subroutine c1f4kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F4KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,4)
  real ( kind = 8 ) ch(in2,l1,4,ido)
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti1
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) tr1
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) wa(ido,3,2)

      if ( 1 < ido ) go to 102
      sn = 1.0D+00 / real ( 4 * l1, kind = 8 )
      if (na == 1) go to 106

      do k = 1, l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         cc(1,k,1,1) = sn*(tr2+tr3)
         cc(1,k,1,3) = sn*(tr2-tr3)
         cc(2,k,1,1) = sn*(ti2+ti3)
         cc(2,k,1,3) = sn*(ti2-ti3)
         cc(1,k,1,2) = sn*(tr1+tr4)
         cc(1,k,1,4) = sn*(tr1-tr4)
         cc(2,k,1,2) = sn*(ti1+ti4)
         cc(2,k,1,4) = sn*(ti1-ti4)
      end do

      return

  106 do 107 k = 1, l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         ch(1,k,1,1) = sn*(tr2+tr3)
         ch(1,k,3,1) = sn*(tr2-tr3)
         ch(2,k,1,1) = sn*(ti2+ti3)
         ch(2,k,3,1) = sn*(ti2-ti3)
         ch(1,k,2,1) = sn*(tr1+tr4)
         ch(1,k,4,1) = sn*(tr1-tr4)
         ch(2,k,2,1) = sn*(ti1+ti4)
         ch(2,k,4,1) = sn*(ti1-ti4)
  107 continue

      return

  102 do 103 k = 1, l1
         ti1 = cc(2,k,1,1)-cc(2,k,1,3)
         ti2 = cc(2,k,1,1)+cc(2,k,1,3)
         tr4 = cc(2,k,1,2)-cc(2,k,1,4)
         ti3 = cc(2,k,1,2)+cc(2,k,1,4)
         tr1 = cc(1,k,1,1)-cc(1,k,1,3)
         tr2 = cc(1,k,1,1)+cc(1,k,1,3)
         ti4 = cc(1,k,1,4)-cc(1,k,1,2)
         tr3 = cc(1,k,1,2)+cc(1,k,1,4)
         ch(1,k,1,1) = tr2+tr3
         ch(1,k,3,1) = tr2-tr3
         ch(2,k,1,1) = ti2+ti3
         ch(2,k,3,1) = ti2-ti3
         ch(1,k,2,1) = tr1+tr4
         ch(1,k,4,1) = tr1-tr4
         ch(2,k,2,1) = ti1+ti4
         ch(2,k,4,1) = ti1-ti4
  103 continue
      do 105 i = 2, ido
         do 104 k = 1, l1
            ti1 = cc(2,k,i,1)-cc(2,k,i,3)
            ti2 = cc(2,k,i,1)+cc(2,k,i,3)
            ti3 = cc(2,k,i,2)+cc(2,k,i,4)
            tr4 = cc(2,k,i,2)-cc(2,k,i,4)
            tr1 = cc(1,k,i,1)-cc(1,k,i,3)
            tr2 = cc(1,k,i,1)+cc(1,k,i,3)
            ti4 = cc(1,k,i,4)-cc(1,k,i,2)
            tr3 = cc(1,k,i,2)+cc(1,k,i,4)
            ch(1,k,1,i) = tr2+tr3
            cr3 = tr2-tr3
            ch(2,k,1,i) = ti2+ti3
            ci3 = ti2-ti3
            cr2 = tr1+tr4
            cr4 = tr1-tr4
            ci2 = ti1+ti4
            ci4 = ti1-ti4
            ch(1,k,2,i) = wa(i,1,1)*cr2+wa(i,1,2)*ci2
            ch(2,k,2,i) = wa(i,1,1)*ci2-wa(i,1,2)*cr2
            ch(1,k,3,i) = wa(i,2,1)*cr3+wa(i,2,2)*ci3
            ch(2,k,3,i) = wa(i,2,1)*ci3-wa(i,2,2)*cr3
            ch(1,k,4,i) = wa(i,3,1)*cr4+wa(i,3,2)*ci4
            ch(2,k,4,i) = wa(i,3,1)*ci4-wa(i,3,2)*cr4
  104    continue
  105 continue

  return
end



subroutine c1f5kf ( ido, l1, na, cc, in1, ch, in2, wa )

!*****************************************************************************80
!
!! C1F5KF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) l1

  real ( kind = 8 ) cc(in1,l1,ido,5)
  real ( kind = 8 ) ch(in2,l1,5,ido)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  real ( kind = 8 ) ci2
  real ( kind = 8 ) ci3
  real ( kind = 8 ) ci4
  real ( kind = 8 ) ci5
  real ( kind = 8 ) cr2
  real ( kind = 8 ) cr3
  real ( kind = 8 ) cr4
  real ( kind = 8 ) cr5
  real ( kind = 8 ) di2
  real ( kind = 8 ) di3
  real ( kind = 8 ) di4
  real ( kind = 8 ) di5
  real ( kind = 8 ) dr2
  real ( kind = 8 ) dr3
  real ( kind = 8 ) dr4
  real ( kind = 8 ) dr5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) ti2
  real ( kind = 8 ) ti3
  real ( kind = 8 ) ti4
  real ( kind = 8 ) ti5
  real ( kind = 8 ), parameter :: ti11 = -0.9510565162951536D+00
  real ( kind = 8 ), parameter :: ti12 = -0.5877852522924731D+00
  real ( kind = 8 ) tr2
  real ( kind = 8 ) tr3
  real ( kind = 8 ) tr4
  real ( kind = 8 ) tr5
  real ( kind = 8 ), parameter :: tr11 =  0.3090169943749474D+00
  real ( kind = 8 ), parameter :: tr12 = -0.8090169943749474D+00
  real ( kind = 8 ) wa(ido,4,2)

      if ( 1 < ido ) go to 102
      sn = 1.0D+00 / real ( 5 * l1, kind = 8 )
      if (na == 1) go to 106

      do k = 1, l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         chold1 = sn*(cc(1,k,1,1)+tr2+tr3)
         chold2 = sn*(cc(2,k,1,1)+ti2+ti3)
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cc(1,k,1,1) = chold1
         cc(2,k,1,1) = chold2
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         cc(1,k,1,2) = sn*(cr2-ci5)
         cc(1,k,1,5) = sn*(cr2+ci5)
         cc(2,k,1,2) = sn*(ci2+cr5)
         cc(2,k,1,3) = sn*(ci3+cr4)
         cc(1,k,1,3) = sn*(cr3-ci4)
         cc(1,k,1,4) = sn*(cr3+ci4)
         cc(2,k,1,4) = sn*(ci3-cr4)
         cc(2,k,1,5) = sn*(ci2-cr5)
      end do

      return

  106 do 107 k = 1, l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         ch(1,k,1,1) = sn*(cc(1,k,1,1)+tr2+tr3)
         ch(2,k,1,1) = sn*(cc(2,k,1,1)+ti2+ti3)
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2,1) = sn*(cr2-ci5)
         ch(1,k,5,1) = sn*(cr2+ci5)
         ch(2,k,2,1) = sn*(ci2+cr5)
         ch(2,k,3,1) = sn*(ci3+cr4)
         ch(1,k,3,1) = sn*(cr3-ci4)
         ch(1,k,4,1) = sn*(cr3+ci4)
         ch(2,k,4,1) = sn*(ci3-cr4)
         ch(2,k,5,1) = sn*(ci2-cr5)
  107 continue

      return

  102 do 103 k = 1, l1
         ti5 = cc(2,k,1,2)-cc(2,k,1,5)
         ti2 = cc(2,k,1,2)+cc(2,k,1,5)
         ti4 = cc(2,k,1,3)-cc(2,k,1,4)
         ti3 = cc(2,k,1,3)+cc(2,k,1,4)
         tr5 = cc(1,k,1,2)-cc(1,k,1,5)
         tr2 = cc(1,k,1,2)+cc(1,k,1,5)
         tr4 = cc(1,k,1,3)-cc(1,k,1,4)
         tr3 = cc(1,k,1,3)+cc(1,k,1,4)
         ch(1,k,1,1) = cc(1,k,1,1)+tr2+tr3
         ch(2,k,1,1) = cc(2,k,1,1)+ti2+ti3
         cr2 = cc(1,k,1,1)+tr11*tr2+tr12*tr3
         ci2 = cc(2,k,1,1)+tr11*ti2+tr12*ti3
         cr3 = cc(1,k,1,1)+tr12*tr2+tr11*tr3
         ci3 = cc(2,k,1,1)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         ch(1,k,2,1) = cr2-ci5
         ch(1,k,5,1) = cr2+ci5
         ch(2,k,2,1) = ci2+cr5
         ch(2,k,3,1) = ci3+cr4
         ch(1,k,3,1) = cr3-ci4
         ch(1,k,4,1) = cr3+ci4
         ch(2,k,4,1) = ci3-cr4
         ch(2,k,5,1) = ci2-cr5
  103 continue

      do 105 i = 2, ido
         do 104 k = 1, l1
            ti5 = cc(2,k,i,2)-cc(2,k,i,5)
            ti2 = cc(2,k,i,2)+cc(2,k,i,5)
            ti4 = cc(2,k,i,3)-cc(2,k,i,4)
            ti3 = cc(2,k,i,3)+cc(2,k,i,4)
            tr5 = cc(1,k,i,2)-cc(1,k,i,5)
            tr2 = cc(1,k,i,2)+cc(1,k,i,5)
            tr4 = cc(1,k,i,3)-cc(1,k,i,4)
            tr3 = cc(1,k,i,3)+cc(1,k,i,4)
            ch(1,k,1,i) = cc(1,k,i,1)+tr2+tr3
            ch(2,k,1,i) = cc(2,k,i,1)+ti2+ti3
            cr2 = cc(1,k,i,1)+tr11*tr2+tr12*tr3
            ci2 = cc(2,k,i,1)+tr11*ti2+tr12*ti3
            cr3 = cc(1,k,i,1)+tr12*tr2+tr11*tr3
            ci3 = cc(2,k,i,1)+tr12*ti2+tr11*ti3
            cr5 = ti11*tr5+ti12*tr4
            ci5 = ti11*ti5+ti12*ti4
            cr4 = ti12*tr5-ti11*tr4
            ci4 = ti12*ti5-ti11*ti4
            dr3 = cr3-ci4
            dr4 = cr3+ci4
            di3 = ci3+cr4
            di4 = ci3-cr4
            dr5 = cr2+ci5
            dr2 = cr2-ci5
            di5 = ci2-cr5
            di2 = ci2+cr5
            ch(1,k,2,i) = wa(i,1,1)*dr2+wa(i,1,2)*di2
            ch(2,k,2,i) = wa(i,1,1)*di2-wa(i,1,2)*dr2
            ch(1,k,3,i) = wa(i,2,1)*dr3+wa(i,2,2)*di3
            ch(2,k,3,i) = wa(i,2,1)*di3-wa(i,2,2)*dr3
            ch(1,k,4,i) = wa(i,3,1)*dr4+wa(i,3,2)*di4
            ch(2,k,4,i) = wa(i,3,1)*di4-wa(i,3,2)*dr4
            ch(1,k,5,i) = wa(i,4,1)*dr5+wa(i,4,2)*di5
            ch(2,k,5,i) = wa(i,4,1)*di5-wa(i,4,2)*dr5
  104    continue
  105 continue

  return
end



subroutine c1fgkf ( ido, ip, l1, lid, na, cc, cc1, in1, ch, ch1, in2, wa )

!*****************************************************************************80
!
!! C1FGKF is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) in1
  integer ( kind = 4 ) in2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) lid

  real ( kind = 8 ) cc(in1,l1,ip,ido)
  real ( kind = 8 ) cc1(in1,lid,ip)
  real ( kind = 8 ) ch(in2,l1,ido,ip)
  real ( kind = 8 ) ch1(in2,lid,ip)
  real ( kind = 8 ) chold1
  real ( kind = 8 ) chold2
  integer ( kind = 4 ) i
  integer ( kind = 4 ) idlj
  integer ( kind = 4 ) ipp2
  integer ( kind = 4 ) ipph
  integer ( kind = 4 ) j
  integer ( kind = 4 ) jc
  integer ( kind = 4 ) k
  integer ( kind = 4 ) ki
  integer ( kind = 4 ) l
  integer ( kind = 4 ) lc
  integer ( kind = 4 ) na
  real ( kind = 8 ) sn
  real ( kind = 8 ) wa(ido,ip-1,2)
  real ( kind = 8 ) wai
  real ( kind = 8 ) war

      ipp2 = ip+2
      ipph = (ip+1)/2
      do ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
      end do

      do 111 j=2,ipph
         jc = ipp2-j
         do 112 ki=1,lid
            ch1(1,ki,j) =  cc1(1,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)-cc1(1,ki,jc)
            ch1(2,ki,j) =  cc1(2,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(2,ki,jc)
  112    continue
  111 continue

      do 118 j=2,ipph
         do 117 ki=1,lid
            cc1(1,ki,1) = cc1(1,ki,1)+ch1(1,ki,j)
            cc1(2,ki,1) = cc1(2,ki,1)+ch1(2,ki,j)
  117    continue
  118 continue

      do 116 l=2,ipph
         lc = ipp2-l
         do 113 ki=1,lid
            cc1(1,ki,l) = ch1(1,ki,1)+wa(1,l-1,1)*ch1(1,ki,2)
            cc1(1,ki,lc) = -wa(1,l-1,2)*ch1(1,ki,ip)
            cc1(2,ki,l) = ch1(2,ki,1)+wa(1,l-1,1)*ch1(2,ki,2)
            cc1(2,ki,lc) = -wa(1,l-1,2)*ch1(2,ki,ip)
  113    continue
         do 115 j=3,ipph
            jc = ipp2-j
            idlj = mod((l-1)*(j-1),ip)
            war = wa(1,idlj,1)
            wai = -wa(1,idlj,2)
            do 114 ki=1,lid
               cc1(1,ki,l) = cc1(1,ki,l)+war*ch1(1,ki,j)
               cc1(1,ki,lc) = cc1(1,ki,lc)+wai*ch1(1,ki,jc)
               cc1(2,ki,l) = cc1(2,ki,l)+war*ch1(2,ki,j)
               cc1(2,ki,lc) = cc1(2,ki,lc)+wai*ch1(2,ki,jc)
  114       continue
  115    continue
  116 continue

      if ( 1 < ido ) go to 136
      sn = 1.0D+00 / real ( ip * l1, kind = 8 )
      if (na == 1) go to 146
      do 149 ki=1,lid
         cc1(1,ki,1) = sn*cc1(1,ki,1)
         cc1(2,ki,1) = sn*cc1(2,ki,1)
  149 continue
      do 120 j=2,ipph
         jc = ipp2-j
         do 119 ki=1,lid
            chold1 = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
            chold2 = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
            cc1(1,ki,j) = chold1
            cc1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
            cc1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
            cc1(1,ki,jc) = chold2
  119    continue
  120 continue
      return

  146 do 147 ki=1,lid
         ch1(1,ki,1) = sn*cc1(1,ki,1)
         ch1(2,ki,1) = sn*cc1(2,ki,1)
  147 continue
      do 145 j=2,ipph
         jc = ipp2-j
         do 144 ki=1,lid
            ch1(1,ki,j) = sn*(cc1(1,ki,j)-cc1(2,ki,jc))
            ch1(2,ki,j) = sn*(cc1(2,ki,j)+cc1(1,ki,jc))
            ch1(1,ki,jc) = sn*(cc1(1,ki,j)+cc1(2,ki,jc))
            ch1(2,ki,jc) = sn*(cc1(2,ki,j)-cc1(1,ki,jc))
  144    continue
  145 continue
      return

  136 do 137 ki=1,lid
         ch1(1,ki,1) = cc1(1,ki,1)
         ch1(2,ki,1) = cc1(2,ki,1)
  137 continue
      do 135 j=2,ipph
         jc = ipp2-j
         do 134 ki=1,lid
            ch1(1,ki,j) = cc1(1,ki,j)-cc1(2,ki,jc)
            ch1(2,ki,j) = cc1(2,ki,j)+cc1(1,ki,jc)
            ch1(1,ki,jc) = cc1(1,ki,j)+cc1(2,ki,jc)
            ch1(2,ki,jc) = cc1(2,ki,j)-cc1(1,ki,jc)
  134    continue
  135 continue
      do 131 i=1,ido
         do 130 k = 1, l1
            cc(1,k,1,i) = ch(1,k,i,1)
            cc(2,k,1,i) = ch(2,k,i,1)
  130    continue
  131 continue
      do 123 j=2,ip
         do 122 k = 1, l1
            cc(1,k,j,1) = ch(1,k,1,j)
            cc(2,k,j,1) = ch(2,k,1,j)
  122    continue
  123 continue
      do 126 j=2,ip
         do 125 i = 2, ido
            do 124 k = 1, l1
               cc(1,k,j,i) = wa(i,j-1,1)*ch(1,k,i,j) &
                            +wa(i,j-1,2)*ch(2,k,i,j)
               cc(2,k,j,i) = wa(i,j-1,1)*ch(2,k,i,j) &
                            -wa(i,j-1,2)*ch(1,k,i,j)
  124       continue
  125    continue
  126 continue

  return
end 



subroutine xerfft ( srname, info )

!*****************************************************************************80
!
!! XERFFT is an error handler for the FFTPACK routines.
!
!  Discussion:
!
!    XERFFT is an error handler for FFTPACK version 5.1 routines.
!    It is called by an FFTPACK 5.1 routine if an input parameter has an
!    invalid value.  A message is printed and execution stops.
!
!    Installers may consider modifying the stop statement in order to
!    call system-specific exception-handling facilities.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, character ( len = * ) SRNAME, the name of the calling routine.
!
!    Input, integer ( kind = 4 ) INFO, an error code.  When a single invalid 
!    parameter in the parameter list of the calling routine has been detected, 
!    INFO is the position of that parameter.  In the case when an illegal 
!    combination of LOT, JUMP, N, and INC has been detected, the calling 
!    subprogram calls XERFFT with INFO = -1.
!
  implicit none

  integer ( kind = 4 ) info
  character ( len = * ) srname

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'XERFFT - Fatal error!'

  if ( 1 <= info ) then
    write ( *, '(a,a,a,i3,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter number ', info, ' had an illegal value.'
  else if ( info == -1 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameters LOT, JUMP, N and INC are inconsistent.'
  else if ( info == -2 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter L is greater than LDIM.'
  else if ( info == -3 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter M is greater than MDIM.'
  else if ( info == -5 ) then
    write ( *, '(a,a,a,a)' ) '  Within ', trim ( srname ), &
      ' input error returned by lower level routine.'
  else if ( info == -6 ) then
    write ( *, '(a,a,a,a)' ) '  On entry to ', trim ( srname ), &
      ' parameter LDIM is less than 2*(L/2+1).'
  end if

  stop
end




subroutine r8_tables ( ido, ip, wa )

!*****************************************************************************80
!
!! R8_TABLES computes trigonometric tables, real double precision arithmetic.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip

  real ( kind = 8 ) arg1
  real ( kind = 8 ) arg2
  real ( kind = 8 ) arg3
  real ( kind = 8 ) arg4
  real ( kind = 8 ) argz
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  real ( kind = 8 ) tpi
  real ( kind = 8 ) wa(ido,ip-1,2)

  tpi = 8.0D+00 * atan ( 1.0D+00 )
  argz = tpi / real ( ip, kind = 8 )
  arg1 = tpi / real ( ido * ip, kind = 8 )

  do j = 2, ip

    arg2 = real ( j - 1, kind = 8 ) * arg1

    do i = 1, ido
      arg3 = real ( i - 1, kind = 8 ) * arg2
      wa(i,j-1,1) = cos ( arg3 )
      wa(i,j-1,2) = sin ( arg3 )
    end do

    if ( 5 < ip ) then
      arg4 = real ( j - 1, kind = 8 ) * argz
      wa(1,j-1,1) = cos ( arg4 )
      wa(1,j-1,2) = sin ( arg4 )
    end if

  end do

  return
end


subroutine r8_factor ( n, nf, fac )

!*****************************************************************************80
!
!! R8_FACTOR factors of an integer for real double precision computations.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the number for which factorization and 
!    other information is needed.
!
!    Output, integer ( kind = 4 ) NF, the number of factors.
!
!    Output, real ( kind = 8 ) FAC(*), a list of factors of N.
!
  implicit none

  real ( kind = 8 ) fac(*)
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  integer ( kind = 4 ) nl
  integer ( kind = 4 ) nq
  integer ( kind = 4 ) nr
  integer ( kind = 4 ) ntry

  nl = n
  nf = 0
  j = 0

  do while ( 1 < nl )

    j = j + 1

    if ( j == 1 ) then
      ntry = 4
    else if ( j == 2 ) then
      ntry = 2
    else if ( j == 3 ) then
      ntry = 3
    else if ( j == 4 ) then
      ntry = 5
    else
      ntry = ntry + 2
    end if

    do

      nq = nl / ntry
      nr = nl - ntry * nq

      if ( nr /= 0 ) then
        exit
      end if

      nf = nf + 1
      fac(nf) = real ( ntry, kind = 8 )
      nl = nq

    end do

  end do

  return
end



subroutine r8_mcfti1 ( n, wa, fnf, fac )

!*****************************************************************************80
!
!! R8_MCFTI1 sets up factors and tables, real double precision arithmetic.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) n
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)
!
!  Get the factorization of N.
!
  call r8_factor ( n, nf, fac )
  fnf = real ( nf, kind = 8 )
  iw = 1
  l1 = 1
!
!  Set up the trigonometric tables.
!
  do k1 = 1, nf
    ip = int ( fac(k1) )
    l2 = l1 * ip
    ido = n / l2
    call r8_tables ( ido, ip, wa(iw) )
    iw = iw + ( ip - 1 ) * ( ido + ido )
    l1 = l2
  end do

  return
end




subroutine cfft1i ( n, wsave, lensav, ier )

!*****************************************************************************80
!
!! CFFT1I: initialization for CFFT1B and CFFT1F.
!
!  Discussion:
!
!    CFFT1I initializes array WSAVE for use in its companion routines 
!    CFFT1B and CFFT1F.  Routine CFFT1I must be called before the first 
!    call to CFFT1B or CFFT1F, and after whenever the value of integer 
!    N changes.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product 
!    of small primes.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Output, real ( kind = 8 ) WSAVE(LENSAV), containing the prime factors 
!    of N and  also containing certain trigonometric values which will be used 
!    in routines CFFT1B or CFFT1F.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    2, input parameter LENSAV not big enough.

  implicit none

  integer ( kind = 4 ) lensav

  integer ( kind = 4 ) ier
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lensav < 2*n + int(log( real ( n, kind = 8 ) )/log( 2.0D+00 )) + 4) then
    ier = 2
    call xerfft ('cfftmi ', 3)
  end if

  if ( n == 1 ) then
    return
  end if

  iw1 = n+n+1

  call r8_mcfti1 ( n, wsave, wsave(iw1), wsave(iw1+1) )

  return
end

subroutine cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

!*****************************************************************************80
!
!! CFFT1F: complex double precision forward fast Fourier transform, 1D.
!
!  Discussion:
!
!    CFFT1F computes the one-dimensional Fourier transform of a single 
!    periodic sequence within a complex array.  This transform is referred 
!    to as the forward transform or Fourier analysis, transforming the 
!    sequence from physical to spectral space.
!
!    This transform is normalized since a call to CFFT1F followed
!    by a call to CFFT1B (or vice-versa) reproduces the original
!    array within roundoff error.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the length of the sequence to be 
!    transformed.  The transform is most efficient when N is a product of 
!    small primes.
!
!    Input, integer ( kind = 4 ) INC, the increment between the locations, in 
!    array C, of two consecutive elements within the sequence to be transformed.
!
!    Input/output, complex ( kind = 8 ) C(LENC) containing the sequence to 
!    be transformed.
!
!    Input, integer ( kind = 4 ) LENC, the dimension of the C array.  
!    LENC must be at least INC*(N-1) + 1.
!
!    Input, real ( kind = 8 ) WSAVE(LENSAV).  WSAVE's contents must be 
!    initialized with a call to CFFT1I before the first call to routine CFFT1F 
!    or CFFT1B for a given transform length N.  WSAVE's contents may be re-used
!    for subsequent calls to CFFT1F and CFFT1B with the same N.
!
!    Input, integer ( kind = 4 ) LENSAV, the dimension of the WSAVE array.  
!    LENSAV must be at least 2*N + INT(LOG(REAL(N))) + 4.
!
!    Workspace, real ( kind = 8 ) WORK(LENWRK).
!
!    Input, integer ( kind = 4 ) LENWRK, the dimension of the WORK array.  
!    LENWRK must be at least 2*N.
!
!    Output, integer ( kind = 4 ) IER, error flag.
!    0, successful exit;
!    1, input parameter LENC   not big enough;
!    2, input parameter LENSAV not big enough;
!    3, input parameter LENWRK not big enough;
!    20, input error returned by lower level routine.
!
  implicit none

  integer ( kind = 4 ) lenc
  integer ( kind = 4 ) lensav
  integer ( kind = 4 ) lenwrk

  complex ( kind = 8 ) c(lenc)
  integer ( kind = 4 ) ier
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) iw1
  integer ( kind = 4 ) n
  real ( kind = 8 ) work(lenwrk)
  real ( kind = 8 ) wsave(lensav)

  ier = 0

  if (lenc < inc*(n-1) + 1) then
    ier = 1
    call xerfft ('cfft1f ', 4)
  else if (lensav < 2*n + int(log( real ( n, kind = 8 ) ) &
    /log( 2.0D+00 )) + 4) then
    ier = 2
    call xerfft ('cfft1f ', 6)
  else if (lenwrk < 2*n) then
    ier = 3
    call xerfft ('cfft1f ', 8)
  end if

  if (n == 1) then
    return
  end if

  iw1 = n + n + 1
  call c1fm1f ( n, inc, c, work, wsave, wsave(iw1), wsave(iw1+1) )

  return
end


subroutine c1fm1f ( n, inc, c, ch, wa, fnf, fac )

!*****************************************************************************80
!
!! C1FM1F is an FFTPACK5.1 auxiliary routine.
!
!  License:
!
!    Licensed under the GNU General Public License (GPL).
!    Copyright (C) 1995-2004, Scientific Computing Division,
!    University Corporation for Atmospheric Research
!
!  Modified:
!
!    15 November 2011
!
!  Author:
!
!    Original FORTRAN77 version by Paul Swarztrauber, Richard Valent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Paul Swarztrauber,
!    Vectorizing the Fast Fourier Transforms,
!    in Parallel Computations,
!    edited by G. Rodrigue,
!    Academic Press, 1982.
!
!    Paul Swarztrauber,
!    Fast Fourier Transform Algorithms for Vector Computers,
!    Parallel Computing, pages 45-63, 1984.
!
!  Parameters:
!
  implicit none

  complex ( kind = 8 ) c(*)
  real ( kind = 8 ) ch(*)
  real ( kind = 8 ) fac(*)
  real ( kind = 8 ) fnf
  integer ( kind = 4 ) ido
  integer ( kind = 4 ) inc
  integer ( kind = 4 ) inc2
  integer ( kind = 4 ) ip
  integer ( kind = 4 ) iw
  integer ( kind = 4 ) k1
  integer ( kind = 4 ) l1
  integer ( kind = 4 ) l2
  integer ( kind = 4 ) lid
  integer ( kind = 4 ) n
  integer ( kind = 4 ) na
  integer ( kind = 4 ) nbr
  integer ( kind = 4 ) nf
  real ( kind = 8 ) wa(*)

      inc2 = inc+inc
      nf = fnf
      na = 0
      l1 = 1
      iw = 1

      do k1=1,nf

         ip = fac(k1)
         l2 = ip*l1
         ido = n/l2
         lid = l1*ido
         nbr = 1+na+2*min(ip-2,4)
         go to (52,62,53,63,54,64,55,65,56,66),nbr
   52    call c1f2kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   62    call c1f2kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   53    call c1f3kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   63    call c1f3kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   54    call c1f4kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   64    call c1f4kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   55    call c1f5kf (ido,l1,na,c,inc2,ch,2,wa(iw))
         go to 120
   65    call c1f5kf (ido,l1,na,ch,2,c,inc2,wa(iw))
         go to 120
   56    call c1fgkf (ido,ip,l1,lid,na,c,c,inc2,ch,ch,2,wa(iw))
         go to 120
   66    call c1fgkf (ido,ip,l1,lid,na,ch,ch,2,c,c,inc2,wa(iw))
  120    l1 = l2
         iw = iw+(ip-1)*(ido+ido)
         if(ip <= 5) na = 1-na
      end do

  return
end

