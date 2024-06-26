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
        rx=rd_x/(2*math.pi)
        ry=rd_y/(2*math.pi)
        if (ax<0)==(ay<0):
          rx=-rx#negative radius means capstan is on the left side of the cable
        else:
          ry=-ry 
        sx.r=rx
        sx.C=ax+1j*math.copysign(rx,ay)
        sy.r=ry 
        sy.C=math.copysign(ry,ax)+1j*ay
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
        acoords = list(zip(*self.anchors))
        self.axes_min = toolhead.Coord(*[min(a) for a in acoords], e=0.)
        self.axes_max = toolhead.Coord(*[max(a) for a in acoords], e=0.)
        self.set_position([0., 0., 0.], ())
        # Setup boundary checks
        max_velocity, max_accel = toolhead.get_max_velocity()
        self.max_z_velocity = config.getfloat('max_z_velocity', max_velocity,
                                              above=0., maxval=max_velocity)
        self.max_z_accel = config.getfloat('max_z_accel', max_accel,
                                           above=0., maxval=max_accel)
        self.steppers=[sx,sy,sz]

    def get_steppers(self):
        return list(self.steppers)
    def calc_position(self, stepper_positions):
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
        Rx=stepper_positions['stepper_x']
        Ry=stepper_positions['stepper_y']
        z=stepper_positions['stepper_z']
        Cx=self.steppers[0].C
        Cy=self.steppers[1].C
        x,y=intersect_circles(Cx.real,0.0,Rx,0.0,Cy.imag,Ry)
        return x,y,z
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
