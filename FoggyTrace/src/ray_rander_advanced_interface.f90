module ray_rander_advanced_interface
    use file_io

    use ray_renderadvanced
    implicit none

    contains

        subroutine test_run()
            implicit none

            integer, parameter :: pixel_x = 640!320
            integer, parameter :: pixel_y = 512!320

            integer, parameter :: no_obj = 100
            integer, parameter :: no_mat = 30

            real,dimension(pixel_x,pixel_y)::img1,img2,img3
            real,dimension(no_obj)::obj
            real,dimension(no_mat)::mat
            integer::anzx,anzy,panz
            real::tof,tofstart,tofend
            real,dimension(pixel_x,pixel_y)::tof1,tof2,tof3

            real::skyint,skycol(3),sunrad,sunint,sunpos(3),suncol(3)

            real,dimension(3)::v1,v2
            real skaletheta

            real :: sphere_centre(3)
            real :: sphere_radius

            real :: mat_rgb(3)
            real :: rs,ws,as,rd,ad
            complex :: refractive_index

            integer :: u
            logical :: rv

            skyint = 0.5 ! [%]
            skycol = (/39.0, 255.0, 210.0/) ! RGB
!            skycol = (/255.0, 0.0, 0.0/) ! RGB

            sunint = 1 ! [%]
            suncol = (/252.0, 255.0, 255.0/) ! RGB
            sunrad = 20 ! [°]
!            sunpos = (/5., 1., 1./) ! [m]
!            sunpos = (/100., 0., 100./) ! [m]
            sunpos = (/2., 0., .5/) ! [m]
            call render_global_par(skyint,skycol,sunint,sunrad,sunpos,suncol)


            v1 = (/1.5, 0., 0.5/)   ! CamFocuspoint
            v2 = (/0., 0., 0.5/)    ! CamFocusLength
            skaletheta = 4!3!1.5
            call raycam_fisheye(v1,v2,skaletheta)

            panz = 100
            anzx = pixel_x
            anzy = pixel_y

            sphere_centre = (/0., -0.3, 0.3/)
            sphere_radius = 0.2
            call rayobject_PutSphere(obj, 1., sphere_centre, sphere_radius)
            sphere_centre = (/-1.0, 0.0, 0.6/)
            call rayobject_PutSphere(obj, 1., sphere_centre, sphere_radius)
            sphere_centre = (/-2.0, 0.3, 0.9/)
            call rayobject_PutSphere(obj, 1., sphere_centre, sphere_radius)
            sphere_centre = (/-3.0, 0.6, 1.2/)
            call rayobject_PutSphere(obj, 1., sphere_centre, sphere_radius)


!            mat_rgb = (/249., 137., 26./)
            mat_rgb = (/255., 255., 255./)
            rs = 0.01    ! Reflektivität specular
            ws = 0.9    ! Wahrscheinlickeit specular
            as = 20     ! Winkel specular
            rd = 0.001    ! Refletivität diffuse
            ad = 50     ! Winkel diffuse
            refractive_index = cmplx(1.33, 0.)
            call raymaterial_PutMaterial(mat, mat_rgb, rs,ws,as,rd,ad, refractive_index)


            call rendermain_advanced(img1,img2,img3, obj,mat, anzx,anzy,panz)
            !call rendermain_advanced(img1,img2,img3, obj,mat, anzx,anzy,panz, tof,tofstart,tofend,tof1,tof2,tof3)

            open(unit=u, file="FoggyTrace_imgRGB.ppm", status="unknown")
            rv = write_ppm_p3(u, img1, img2, img3, rescale=.false.)
            close(u)

            open(unit=u, file="FoggyTrace_imgRGB_scaled.ppm", status="unknown")
            rv = write_ppm_p3(u, img1, img2, img3, rescale=.true.)
            close(u)

            open(unit=u, file="FoggyTrace_imgRGB_hdr.ppm", status="unknown")
            rv = write_ppm_p3(u, img1, img2, img3, hdr=.true.)
            close(u)

        end subroutine
end module ray_rander_advanced_interface
