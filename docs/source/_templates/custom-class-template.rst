{{ fullname | escape | underline}}

- full name: {{ fullname | escape }}
- parent module: :mod:`{{ module }}`
- type: {{ objtype }}

.. currentmodule:: {{ module }}

.. rubric:: Inheritance Diagram

.. inheritance-diagram:: {{ fullname }}
    :parts: 1

|

.. autoclass:: {{ objname }}
   :members:
   :inherited-members:
   :show-inheritance:
   :special-members: __init__

   {% block methods %}

   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
