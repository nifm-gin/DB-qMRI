classdef sigmoidLayer < nnet.layer.Layer
    methods
        function layer = sigmoidLayer(name) 
            % Set layer name
            if nargin == 2
                layer.Name = name;
            end
            % Set layer description
            layer.Description = 'sigmoidLayer'; 
        end
        function Z = predict(layer,X)
            % Forward input data through the layer and output the result
            Z = exp(X)./(exp(X)+1);
        end
        function dLdX = backward(layer, X ,Z,dLdZ,memory)
            % Backward propagate the derivative of the loss function through 
            % the layer 
            dLdX = Z.*(1-Z) .* dLdZ;
        end
    end
 end