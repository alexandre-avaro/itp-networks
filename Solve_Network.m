function [I_Array, V_Array]=Solve_Network(Channel_Array, Node_Array, Num_Edges, Num_Nodes, Network_V, Network_I)
% here we solve for the current and voltage values within the network
% we construct a matrix which contains the equations-->the unknown vector
% is stacked with I values first (for edges) and then V values (for nodes)

fixed_Voltage=false;
fixed_Current=true;


A=zeros(Num_Edges+Num_Nodes,Num_Edges+Num_Nodes);
RHS=zeros(Num_Edges+Num_Nodes,1);

if(fixed_Voltage)
    for i=1:Num_Nodes
        if(i==1)
            A(i,Num_Edges+i)=1;
            RHS(i)=0;
        elseif(i==Num_Nodes)
            A(i,Num_Edges+i)=1;
            RHS(i)=Network_V;
        else
            for j=1:Node_Array(i).Num_Edges_L
                A(i, Node_Array(i).Edges_L(j))=1;
            end
            for j=1:Node_Array(i).Num_Edges_R
                A(i, Node_Array(i).Edges_R(j))=-1;
            end
            RHS(i)=0;
        end
    end
    
    for i=1:Num_Edges
        A(Num_Nodes+i, i)=-Channel_Array(i).R;
        A(Num_Nodes+i, Num_Edges+Channel_Array(i).Node_L)=-1;
        A(Num_Nodes+i, Num_Edges+Channel_Array(i).Node_R)=1;
        RHS(Num_Nodes+i)=0;
    end
end

if(fixed_Current)
    for i=1:Num_Nodes
        if(i==1)
            A(i,Num_Edges+i)=1;
            RHS(i)=0;
        elseif(i==Num_Nodes)
            for j=1:Node_Array(i).Num_Edges_L
                A(i, Node_Array(i).Edges_L(j))=1;
            end
            RHS(i)=Network_I;
        else
            for j=1:Node_Array(i).Num_Edges_L
                A(i, Node_Array(i).Edges_L(j))=1;
            end
            for j=1:Node_Array(i).Num_Edges_R
                A(i, Node_Array(i).Edges_R(j))=-1;
            end
            RHS(i)=0;
        end
    end
    
    for i=1:Num_Edges
        A(Num_Nodes+i, i)=-Channel_Array(i).R;
        A(Num_Nodes+i, Num_Edges+Channel_Array(i).Node_L)=-1;
        A(Num_Nodes+i, Num_Edges+Channel_Array(i).Node_R)=1;
        RHS(Num_Nodes+i)=0;
    end
end
Solution=A\RHS;
I_Array=Solution(1:Num_Edges);
V_Array=Solution(Num_Edges+1:Num_Edges+Num_Nodes);
end


