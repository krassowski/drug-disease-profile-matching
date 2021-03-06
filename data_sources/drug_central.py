from pandas import DataFrame
from sqlalchemy.ext.automap import automap_base
from sqlalchemy import create_engine, MetaData
from sqlalchemy_utils import generic_repr
from sqlalchemy.orm import scoped_session, sessionmaker
from thirdparty.sqlalchemy import name_for_scalar_relationship, pluralize_collection, camelize_classname

m = MetaData()
Base = generic_repr(automap_base(metadata=m))


def _repr_html_(self):
    from helpers.gui.namespace import HorizontalNamespace
    data = {
        k: v
        for k, v in self.__dict__.items()
        if not k.startswith('_')
    }
    namespace = HorizontalNamespace(data)
    namespace.__name__ = self.__class__.__name__
    return namespace._repr_html_()


Base._repr_html_ = _repr_html_

# ./data/drugcentral/download.sh
# ./data/drugcentral/install.sh
engine = create_engine("postgresql+psycopg2://drugcentral:password@localhost/drugcentral")

m.reflect(engine, views=True)

# reflect the tables
omop_relationship_doid_view = m.tables['omop_relationship_doid_view']


class OmopDoid(Base):
    __table__ = omop_relationship_doid_view
    __mapper_args__ = {
        'primary_key': [omop_relationship_doid_view.c.id]
    }


Base.prepare(
    engine,
    name_for_scalar_relationship=name_for_scalar_relationship,
    classname_for_table=camelize_classname,
    name_for_collection_relationship=pluralize_collection
)

session = scoped_session(sessionmaker(engine))
Base.query = session.query_property()


locals().update(Base.classes.items())
Structure = Base.classes.Structures


def summarize_drugs_for_concept(query, relation=None, approvals=False, quantity=False):
    drugs = []
    if isinstance(query, str):
        concept_filter = (OmopRelationship.concept_name == query)
    else:
        concept_filter = query

    for omop in OmopRelationship.query.filter(concept_filter):

        if relation and relation != omop.relationship_name:
            continue

        structure = Structure.query.filter_by(id=omop.struct_id).one()

        for use_in_a_drug in structure.active_ingredients:
            drug = use_in_a_drug.product

            row = {
                'relation': omop.relationship_name,
                'disease': omop.concept_name,
                'structure': structure.name,
                'product_name': drug.product_name,
                'generic_name': drug.generic_name,
                'route': drug.route,
                'orhter_active_ingredients': drug.active_ingredient_count - 1,
                'marketing_status': drug.marketing_status,
            }
            if quantity:
                row.udpate({
                    'quantity': use_in_a_drug.quantity,
                    'unit': use_in_a_drug.unit
                })
            if approvals:
                row['approvals'] = structure.approvals

            drugs.append(row)

    return DataFrame(drugs).drop_duplicates()


def all_names_for(*substances, exclude_given_names=False):
    return {
        synonym.lname
        for substance_name in substances
        for original_onym in Synonyms.query.filter_by(name=substance_name)
        for synonym in Synonyms.query.filter_by(id=original_onym.id)
        if not exclude_given_names or synonym.lname not in substances
    }


def substances_for(query, contra=False, verbose=False):
    """
    query: e.g. 'Carcinoma of female breast', use '%' for fuzzy matching
    """
    results = set()
    diseases = set()

    query = OmopDoid.concept_name.ilike(query) if '%' in query else OmopDoid.concept_name == query
    relationship = 'contraindication' if contra else 'indication'

    for rel in OmopDoid.query.filter_by(relationship_name=relationship).filter(query):
        diseases.add(rel.concept_name)
        new_result = Structure.query.filter_by(id=rel.struct_id).one().name
        if verbose:
            if new_result in results:
                print(f'Intersection with previous results {new_result}')
        results.add(new_result)
    if verbose:
        print(diseases)
    return results


__all__ = [
    'Base',
    'session',
    'OmopDoid',
    'Structure',
    'summarize_drugs_for_concept',
    'substances_for',
    *Base.classes.keys()
]
