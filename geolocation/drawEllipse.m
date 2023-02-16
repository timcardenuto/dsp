function eep = drawEllipse(xhat,P)
    % Create and scale ellipse
    [eigenvectors, eigenvalues] = eig(P);
    k = 5.9915;
    theta = linspace(0,2*pi);
    ellipse = [sqrt(k*max(eigenvalues(1),eigenvalues(4)))*cos(theta);...
        sqrt(k*min(eigenvalues(1),eigenvalues(4)))*sin(theta)]';
    % Rotate ellipse
    eigenvector = eigenvectors(:,2);
    angle = atan2(eigenvector(2), eigenvector(1));
    if(angle < 0)
        angle = angle + 2*pi;
    end
    R = [cos(angle),sin(angle);-sin(angle),cos(angle)];
    ellipse = ellipse * R;
    % Shift ellipse and plot
    plot(xhat(1), xhat(2), '.b')
    hold on
    eep = plot((ellipse(:,1)+xhat(1)),(ellipse(:,2)+xhat(2)),'-b')
end