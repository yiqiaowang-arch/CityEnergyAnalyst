from fastapi import APIRouter

import cea.interfaces.dashboard.api.inputs as inputs
import cea.interfaces.dashboard.api.contents as contents
import cea.interfaces.dashboard.api.dashboard as dashboard
import cea.interfaces.dashboard.api.databases as databases
import cea.interfaces.dashboard.api.glossary as glossary
import cea.interfaces.dashboard.api.project as project
import cea.interfaces.dashboard.api.tools as tools

router = APIRouter()

router.include_router(inputs.router, prefix="/inputs")
router.include_router(contents.router, prefix="/contents")
router.include_router(dashboard.router, prefix="/dashboard")
router.include_router(databases.router, prefix="/databases")
router.include_router(glossary.router, prefix="/glossary")
router.include_router(project.router, prefix="/project")
router.include_router(tools.router, prefix="/tools")
