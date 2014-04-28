load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	x = fspan(1, 4, 4)
	errors = (/0.00262,0.00524,0.00717,0.01071/)

	; reference Lines
	ref_line_1st = 0.01*x
	ref_line_2nd = 0.0001*x^2

	wks = gsn_open_wks("x11", "convergence")

	res = True
	res@gsnDraw = False
	res@gsnFrame = False
	res@vpWidthF = 0.4
	res@vpHeightF = 0.3
	res@xyXStyle = "Log"
	res@xyYStyle = "Log"
	res@xyLineThicknessF = 5.0
	res@trXReverse = True
	res@trXMinF = 0.5
	res@trXMaxF = 4.5
	res@trYMinF = 1.0e-8
	res@trYMaxF = 1
	res@tmXBTickStartF = res@trXMinF
	res@tmXBTickEndF = res@trXMaxF
	res@tmYLTickStartF = 1.0e-8
	res@tmYLTickEndF = 1
	res@tmXMajorGrid = True
	res@tmYMajorGrid = True
	res@tmXMajorGridThicknessF = 1.0
	res@tmYMajorGridThicknessF = 1.0
	res@tmXMajorGridLineDashPattern = 2
	res@tmYMajorGridLineDashPattern = 2

	res@xyMarkLineMode = "MarkLines"
	res@xyMarkerSizeF = 10
	plot0 = gsn_csm_xy(wks, x, errors, res)

	res@xyMarkLineMode = "Lines"
	plot1 = gsn_csm_xy(wks, x, ref_line_1st, res)
	plot2 = gsn_csm_xy(wks, x, ref_line_2nd, res)

	overlay(plot0, plot1)
	overlay(plot0, plot2)

	draw(plot0)
	frame(wks)

end