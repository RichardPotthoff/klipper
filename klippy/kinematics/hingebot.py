# Code for handling the kinematics of cable winch robots
#
# Copyright (C) 2018-2021  Kevin O'Connor <kevin@koconnor.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.
import stepper, mathutil, math

class HingebotKinematics:
    def __init__(self, toolhead, config):
        # Setup steppers at each anchor
        stepper_config = config.getsection('stepper_x')
        sx = stepper.PrinterStepper(stepper_config)
        rd_x=stepper_config.getfloat('rotation_distance')
        ax=stepper_config.getfloat('anchor')
        stepper_config = config.getsection('stepper_y')
        sy = stepper.PrinterStepper(stepper_config)        
        rd_y=stepper_config.getfloat('rotation_distance')
        ay=stepper_config.getfloat('anchor')
        stepper_config = config.getsection('stepper_z')
        sz = stepper.PrinterStepper(stepper_config)
        self.steppers=[sx,sy,sz]
        rx=rd_x/(2*math.pi)
        ry=rd_y/(2*math.pi)
        if (ax<0)==(ay<0):
          rx=-rx#negative radius means capstan is on the left side of the cable
        else:
          ry=-ry 
        sx.r=rx # capstan radius 
        sx.C=ax+1j*math.copysign(rx,ay) #center coordinates of capstan
        sx.ang0=math.atan2(-sx.C.imag,-sx.C.real) #angle towards origin 
        sy.r=ry 
        sy.C=math.copysign(ry,ax)+1j*ay
        sy.ang0=math.atan2(-sy.C.imag,-sy.C.real)
        sx.setup_itersolve('hingebot_stepper_alloc',sx.C.real,sx.C.imag,0.0,sx.r)
        sx.set_trapq(toolhead.get_trapq())
        toolhead.register_step_generator(sx.generate_steps)
        sy.setup_itersolve('hingebot_stepper_alloc',sy.C.real,sy.C.imag,0.0,sy.r)
        sy.set_trapq(toolhead.get_trapq())
        toolhead.register_step_generator(sy.generate_steps)
        sz.setup_itersolve('cartesian_stepper_alloc','z'.encode()) 
        sz.set_trapq(toolhead.get_trapq())
        toolhead.register_step_generator(sz.generate_steps)
        # Setup boundary checks
        self.axes_min = toolhead.Coord(-abs(ax),-abs(ay),-1.0, e=0.)
        self.axes_max = toolhead.Coord(abs(ax),abs(ay),1000.0, e=0.)
        self.set_position([0., 0., 0.], ())
        # Setup boundary checks
        max_velocity, max_accel = toolhead.get_max_velocity()
        self.max_z_velocity = config.getfloat('max_z_velocity', max_velocity,
                                              above=0., maxval=max_velocity)
        self.max_z_accel = config.getfloat('max_z_accel', max_accel,
                                           above=0., maxval=max_accel)
        

    def get_steppers(self):
        return list(self.steppers)
    def calc_position(self, stepper_positions,p_guess=0.0+0.0j,max_it=10,eps=1e-6):
        from math import sqrt, pi, sin, cos, atan2, copysign
        def evolvent(t,t0=0.0,r=1.0):
            tt0=t+t0
            return r*(cos(tt0)+t*sin(tt0)+1j*(sin(tt0)-t*cos(tt0)))
        def sss(c,a,b):#triangle with 3 sides: return angle opposite first side 'c'
            cosgamma=(a*a+b*b-c*c)/(2*a*b)
            return cosgamma 
        def intersect_circles(c1x,c1y,r1,c2x,c2y,r2):# intersection point of 2 circles
            dx,dy=c1x-c2x,c1y-c2y
            dl=sqrt(dx*dx+dy*dy)
            ex,ey=dx/dl,dy/dl 
            cos1=sss(r1,dl,r2)
            sin1=sqrt(1.0-cos1*cos1)
            p1x,p1y=c2x+r2*(ex*cos1+ey*sin1),c2y+r2*(-ex*sin1+ey*cos1)
            p2x,p2y=c2x+r2*(ex*cos1-ey*sin1),c2y+r2*(ex*sin1+ey*cos1)
            if (p1x*p1x+p1y*p1y)<(p2x*p2x+p2y*p2y):#point closest to the origin
                return p1x,p1y
            else:
                return p2x,p2y
        spx=stepper_positions['stepper_x']
        spy=stepper_positions['stepper_y']
        z=stepper_positions['stepper_z']
        pe=[None]*2#cable end point
        pt=[None]*2#cable tangent point (at capstan)
        ucl=[None]*2#unwound cable length
        for j in range(max_it):
          for i,(stepper,sp) in enumerate(zip(self.steppers[:2],[spx,spy])):
            drot1=(p_guess-stepper.C)*((abs(p_guess-stepper.C)**2-
                   stepper.r**2)**0.5+1j*stepper.r)
            pt[i] = -1j*stepper.r * drot1/abs(drot1) + stepper.C #tangent point = center of approx. circle
            drot=drot1*(0+0j-stepper.C).conjugate()#rotation relative to stepper
            ucl[i]=sp-stepper.r*atan2(drot.imag,drot.real)# unwound cable length = radius of approx. circle
            pe[i]=evolvent(t=-ucl[i]/stepper.r, 
                    t0=sp/stepper.r+stepper.ang0-copysign(pi/2,stepper.r), 
                    r=abs(stepper.r)
                    ) + stepper.C
          newx,newy=intersect_circles(pt[0].real,pt[0].imag,ucl[0],
                      pt[1].real,pt[1].imag,ucl[1])
          p_guess=newx+1j*newy
          error=abs(pe[0]-pe[1])
          if error<eps: break
        return newx,newy,z
    def set_position(self, newpos, homing_axes):
        for s in self.steppers:
            s.set_position(newpos)
    def home(self, homing_state):
        # XXX - homing not implemented
        homing_state.set_axes([0, 1, 2])
        homing_state.set_homed_position([0., 0., 0.])
    def check_move(self, move):
        # XXX - boundary checks and speed limits not implemented
        #pass
        if not move.axes_d[2]:
            # Normal XY move - use defaults
            return
        z_ratio = move.move_d / abs(move.axes_d[2])
        move.limit_speed(
            self.max_z_velocity * z_ratio, self.max_z_accel * z_ratio)
    def get_status(self, eventtime):
        # XXX - homed_checks and rail limits not implemented
        return {
            'homed_axes': 'xyz',
            'axis_minimum': self.axes_min,
            'axis_maximum': self.axes_max,
        }

def load_kinematics(toolhead, config):
    return HingebotKinematics(toolhead, config)
