class Particle:
    
    def __init__(self, mass, x, y, vx, vy):
        self.mass = mass
        self.x = x
        self.y = y
        self.vx = vx
        self.vy = vy
        self.fx = 0
        self.fy = 0
        self.G = 6.673e-11
    
    ##resets the net force vectors on the given mass after each time step
    def resetForce(self):
        self.fx = 0
        self.fy = 0
    
    ##updates the attributes of the given mass object after each time step
    def update(self, dt):
        self.vx += dt * self.fx/self.mass
        self.vy += dt * self.fy/self.mass
        self.x += dt * self.vx
        self.y += dt * self.vy
        
    ##gets the force of the parent star at this point in the orbit 
    def getForce(self, particle):
        dx = self.x - particle.x
        dy = self.y - particle.y
        d = ((dx ** 2) + (dy ** 2)) ** 0.50
        Fg =  -self.G * self.mass * particle.mass/((d+3e4)**2)
        self.fx += Fg * dx/d
        self.fy += Fg * dy/d
        
    def getAvgVel(self, star):
        radius = (((self.x-star.x)**2) + ((self.y-star.y)**2))** 0.50
        PI = 3.14159265358979
        period = ((4*(PI ** 2)*(radius **3))/(self.G * star.mass)) ** 0.50
        circumference = 2 * PI * radius
        return circumference/period
    
    
    
        