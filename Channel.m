classdef Channel < handle % handle class so properties persist
    properties
        Length
        Area
        R
        Sigma_T=0.5
        Sigma_L=0.9
        Sigma_i=(0.9+0.5)/2
        Mu_T=1.95e-9
        Mu_L=79.1e-9
        Volume=0
        Peak_x=[0]
        Peak_V=[0]
        Has_Peak=false
        Num_Peaks=0
        Node_L
        Node_R
        X_L=0
        Y_L=0
        X_R=0
        Y_R=0
    end

    methods
        function obj = Channel(Length, Area, Node_L, Node_R, Has_Peak, Num_Peaks, Peak_x)
            obj.Length = Length;
            obj.Area = Area;
            obj.Node_L = Node_L;
            obj.Node_R = Node_R;
            obj.Has_Peak = Has_Peak;
            obj.Num_Peaks = Num_Peaks;
            obj.Peak_x = Peak_x;
        end

        function Compute_Volume(obj)
            obj.Volume = obj.Length * obj.Area;
        end

        function Compute_R(obj)
            if obj.Num_Peaks < 2
                obj.R = obj.Peak_x(1)/obj.Area*(1/obj.Sigma_T-1/obj.Sigma_L)+obj.Length/(obj.Sigma_L*obj.Area);
            else
                obj.R = obj.Peak_x(2)/obj.Area*(1/obj.Sigma_T-1/obj.Sigma_i)+obj.Peak_x(1)/obj.Area*(1/obj.Sigma_i-1/obj.Sigma_L)+obj.Length/(obj.Sigma_L*obj.Area);
            end 
        end
    end
end