function [central_model] = build_central(select, A, lhs, rhs, ub, lb, f)
%SOLVE_CENTRAL_AND_WRITE_VARIABLES

central_model = struct;
central_model.Model.A = A;
central_model.Model.lhs = lhs;
central_model.Model.rhs = rhs;
central_model.Model.ub = ub;
central_model.Model.lb = lb;
central_model.Model.obj = f;
end

