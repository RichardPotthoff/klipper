// Cable winch stepper kinematics
//
// Copyright (C) 2018-2019  Kevin O'Connor <kevin@koconnor.net>
//
// This file may be distributed under the terms of the GNU GPLv3 license.

#include <math.h> // sqrt, atan2
#include <stddef.h> // offsetof
#include <stdlib.h> // malloc
#include <string.h> // memset
#include "compiler.h" // __visible
#include "itersolve.h" // struct stepper_kinematics
#include "trapq.h" // move_get_coord

struct hingebot_stepper {
    struct stepper_kinematics sk;
    struct coord anchor;
    double r;
    double ang0;
};

static double
hingebot_stepper_calc_position(struct stepper_kinematics *sk, struct move *m
                            , double move_time)
{
    struct hingebot_stepper *hs = container_of(sk, struct hingebot_stepper, sk);
    struct coord c = move_get_coord(m, move_time);
    double dx = hs->anchor.x - c.x, dy = hs->anchor.y - c.y;
    double ucl = sqrt(dx*dx+dy*dy-hs->r*hs->r); // length of unwound cable
    double dl_ang=atan2(-dy,-dx) - hs->ang0; //angle difference between capstan->origin and capstan->point
    if (dl_ang < -M_PI) {dl_ang+=2*M_PI;}
    else if (dl_ang > M_PI) {dl_ang-=2*M_PI;}
    double rl_ang=atan2(hs->r,ucl);
    return hs->r*(dl_ang+rl_ang) + ucl;
}

struct stepper_kinematics * __visible
hingebot_stepper_alloc(double anchor_x, double anchor_y, double anchor_z, double r)
{
    struct hingebot_stepper *hs = malloc(sizeof(*hs));
    memset(hs, 0, sizeof(*hs));
    hs->anchor.x = anchor_x;
    hs->anchor.y = anchor_y;
    hs->anchor.z = anchor_z;
    hs->r=r;
    hs->ang0=atan2(-anchor_y,-anchor_x);
    hs->sk.calc_position_cb = hingebot_stepper_calc_position;
    hs->sk.active_flags = AF_X | AF_Y ;//| AF_Z;
    return &hs->sk;
}
