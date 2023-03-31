ha = phased.URA(10);
ha.Element.BackBaffled = true;
 plotResponse(ha,60e6,60e6,'RespCut','3D','Format','Polar')
