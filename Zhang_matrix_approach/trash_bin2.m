%trash2            
%             fdvx_dksi=WFdksi(i,j,vx);
%             fdvy_dksi=WFdksi(i,j,vy); 
%             fsigmaxx_dksi=WFdksi(i,j,sigmaxx);
%             fsigmayy_dksi=WFdksi(i,j,sigmayy);
%             fsigmaxy_dksi=WFdksi(i,j,sigmaxy);
%             Uf_dksi=[fdvx_dksi, fdvy_dksi, fdsigmaxx_dksi, fdsigmayy_dksi, fdsigmaxy_dksi]';
%                         
%             bdvx_dksi=WBdksi(i,j,vx);
%             bdvy_dksi=WBdksi(i,j,vy);
%             bsigmaxx_dksi=WBdksi(i,j,sigmaxx);
%             bsigmayy_dksi=WBdksi(i,j,sigmayy);
%             bsigmaxy_dksi=WBdksi(i,j,sigmaxy);
%             Ub_dksi=[bdvx_dksi, bdvy_dksi, bdsigmaxx_dksi, bdsigmayy_dksi, bdsigmaxy_dksi]';
%             
%             %-------------------------------------------------
%             %dEta
%             fdvx_deta=WFdeta(i,j,vx);
%             fdvy_deta=WFdeta(i,j,vy);
%             fsigmaxx_deta=WFdeta(i,j,sigmaxx);
%             fsigmayy_deta=WFdeta(i,j,sigmayy);
%             fsigmaxy_deta=WFdeta(i,j,sigmaxy);
%             Uf_deta=[fdvx_deta, fdvy_deta, fdsigmaxx_deta, fdsigmayy_deta, fdsigmaxy_deta]';
%                         
%             bdvx_deta=WBdeta(i,j,vx);
%             bdvy_deta=WBdeta(i,j,vy);
%             bsigmaxx_deta=WBdeta(i,j,sigmaxx);
%             bsigmayy_deta=WBdeta(i,j,sigmayy);
%             bsigmaxy_deta=WBdeta(i,j,sigmaxy);
%             Ub_deta=[bdvx_deta, bdvy_deta, bdsigmaxx_deta, bdsigmayy_deta, bdsigmaxy_deta]';
%             
