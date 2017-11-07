################################################################
#Function to Solve Quadratic Equation (positive root)
###############################################################

quad_roots<-function(a,b,c){
  c(((-b-sqrt(b^2-4*a*c))/(2*a)),((-b+sqrt(b^2-4*a*c))/(2*a)))
}
