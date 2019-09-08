$title

$onText
Problem formulation
$offText
*------------------------------------------------------
*$include "problem10.gms"
*$include "problem20.gms"
*$include "problem30.gms"
*$include "problem40.gms"
*$include "problem50.gms"
*$include "problem60.gms"
*$include "problem70.gms"
*------------------------------------------------------

Alias (n,i,t);

*Scalar BigNumber;

Sets
t 'Periods' /p1*p30/
i 'Suppliers' /s1*s60/
pp 'Price'/c1*c10000/;

Scalar OverallCost 'maximum total cost for dynamic programming';

Parameters
Cost(i,t) 'per unit delivery cost from supplier i in period t'
LB(i,t)   'lower bound demand from supplier i in period t'
UB(i,t)   'lower bound demand from supplier i in period t'
Demand(t) 'summary demand in period t'

Integer Variables
x(i,t),
phi(i,v),
y(t);

Binary Variable z(i,t);

Equations
    mainlotsizing(t)          'ordinary lot-sizing constraint'
    nobacktracking(t)         'no-backtracking constraint'
    lowersupplieropening(i,t) 'open or not supply from supplier i at period t'
    uppersupplieropening(i,t) 'open or not supply from supplier i at period t'
    obj                       'total cost objective';


mainlotsizing(t)$(t>1)..y(t) =e= y(t-1)+sum(i, x(i,t))-Demand(t);
nobacktracking(t)..y(t-1)+sum(i, x(i,t)) =g= Demand(t);
lowersupplieropening(i,t)..x(i,t) =g= LB(i,t)*z(i,t);
uppersupplieropening(i,t)..x(i,t) =l= UB(i,t)*z(i,t);
obj..sum((i,t), Cost(i,t)*x(i,t));

Model lotsizing / all /;

solve lotsizing using mip min obj;
$onText
Решили исходную задачу CPLEX, далее решаем эвристическим алгоритмом:
1. Алгоритм ДП, решающий задачу с краевыми значениями. В случае
точечных интервалов дает оптимум задачи.
2. Начинаем снижаться от этого решения, соотнося этот процесс с
леммой о виде оптимального решения.
$offText
*либо все с нуля, либо с +inf
loop(t,
    loop(i,
        loop(v,x_d(i,t,v)= +inf;
        );
    );
);


*Прямой ход динамического программирования
loop(t,
    loop(i$(v>sum(t, Demand(t))),
        if(i<card(i),
        x_d(i,t,v)=min( x_d(i-1,t,v), Cost(i,t)*LB(i,t)+x_d(i-1,t,v-LB(i,t)), Cost(i,t)*UB(i,t)+x_d(i-1,t,v-UB(i,t)));
        continue;
        );
        x_d(i,t,v)=min( x_d(i-1,t,v), Cost(i,t)*LB(i,t)+x_d(i-1,t,v-LB(i,t)), Cost(i,t)*UB(i,t)+x_d(i-1,t,v-UB(i,t)));
    );
);
*Обратный ход динамического программирования: построение полученного ДП решения

*Покоординатный жадный спуск
*Выбираем координату с ?наименьшим? превышением порога допустимости по суммарному объему
