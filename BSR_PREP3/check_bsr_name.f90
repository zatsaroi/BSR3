!======================================================================
      Subroutine Check_BSR_name(AF)
!======================================================================
!     rename old-version file-name as 'name.c' just to 'name'
!----------------------------------------------------------------------
      Character(*) :: AF
      i = LEN_TRIM(AF)
      if(i.gt.2.and.AF(i-1:i).eq.'.c') AF(i-1:i) = '  '

      End Subroutine Check_BSR_name
