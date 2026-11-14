function md = actual_flow(mcmd, z, p)
%ACTUAL_FLOW  Saturated actuator: proportional loop error clipped to [0, mMax].
md = min(p.mMax, max(0, mcmd.*(1 + p.flowErrRelSd*z)));
end
